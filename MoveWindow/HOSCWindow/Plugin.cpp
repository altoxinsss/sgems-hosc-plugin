
#define NOMINMAX
#include <algorithm>
#include <cmath>
#include <GsTLAppli/geostat/parameters_handler.h>
#include <GsTLAppli/utils/string_manipulation.h>
#include <GsTLAppli/utils/error_messages_handler.h>
#include <GsTLAppli/utils/gstl_plugins.h>
#include <GsTLAppli/geostat/utilities.h>
#include <GsTLAppli/appli/utilities.h>
#include <GsTLAppli/grid/grid_model/grid_initializer.h>
#include <GsTLAppli/grid/grid_model/sgrid_cursor.h>
#include <GsTLAppli/grid/grid_model/rgrid.h>
#include <GsTLAppli/grid/grid_model/point_set.h>
#include <GsTLAppli/grid/grid_model/gval_iterator.h>
#include <GsTL/math/math_functions.h>
#include <GsTLAppli/math/gstlpoint.h>
#include <GsTLAppli/grid/grid_model/cartesian_grid.h>
#include <GsTLAppli/appli/project.h>
#include <GsTLAppli/grid/grid_model/geostat_grid.h>
#include <GsTLAppli/appli/manager_repository.h>

#include <GsTLAppli/actions/common.h>
#include <GsTLAppli/actions/action.h>

#include <ctype.h>
#include <direct.h>
#include <exception>
#include <errno.h>
#include <limits.h>
#include <new>
#include <stdarg.h>
#include <sstream>
#include <fstream>
#include <time.h>
#include <excpt.h>
#include <windows.h>
// Windows headers define `interface` as a macro for COM declarations.  The
// legacy SGeMS API has a real `interface()` member function, so preserve that
// C++ identifier after obtaining the Win32 process/thread APIs.
#ifdef interface
#undef interface
#endif

#include "Plugin.h"
#include "PCA.h"
#include "../../shared/ValidationFeedback.h"

namespace wmc = weighted_moment_core;

namespace {
static volatile LONG hosc_window_trace_sequence = 0;
static volatile LONG hosc_window_run_sequence = 0;

unsigned long next_hosc_window_run_sequence()
{
	return (unsigned long)InterlockedIncrement(&hosc_window_run_sequence);
}

std::string parameter_value_compat(const Parameters_handler* parameters, const char* name)
{
    if(parameters == NULL || name == NULL)
    {
        return std::string();
    }

    // Parameters_handler is the stable host boundary.  Accessing the XML
    // implementation's Qt DOM members directly tied this plugin to a private
    // object layout and required unavailable Qt headers.
    return parameters->value(std::string(name));
}

void append_hosc_window_trace(const char* message)
{
    FILE* trace = fopen("HOSC_WINDOW_debug.log", "a");
    if(trace == NULL)
    {
        return;
    }
	fprintf(
		trace,
		"%ld pid=%lu tid=%lu trace_sequence=%ld %s\n",
		(long)time(NULL),
		(unsigned long)GetCurrentProcessId(),
		(unsigned long)GetCurrentThreadId(),
		(long)InterlockedIncrement(&hosc_window_trace_sequence),
		message == NULL ? "" : message);
    fclose(trace);
}

void append_hosc_window_tracef(const char* format, ...)
{
    char buffer[2048];
    va_list args;
    va_start(args, format);
#if defined(_MSC_VER)
    _vsnprintf(buffer, sizeof(buffer) - 1, format, args);
#else
    vsnprintf(buffer, sizeof(buffer) - 1, format, args);
#endif
    buffer[sizeof(buffer) - 1] = '\0';
    va_end(args);
    append_hosc_window_trace(buffer);
}

void append_hosc_window_instance_trace(
	const void* instance,
	unsigned long run_sequence,
	const char* stage)
{
	append_hosc_window_tracef(
		"run_sequence=%lu instance=%p stage=%s",
		run_sequence,
		instance,
		stage == NULL ? "" : stage);
}

std::string normalize_algorithm_text(const std::string& text)
{
    std::string normalized = text;
    for(size_t i = 0; i < normalized.size(); i++)
    {
        normalized[i] = (char)tolower((unsigned char)normalized[i]);
    }
    return normalized;
}

Algorithm parse_algorithm_selection(const std::string& text)
{
    std::string normalized = normalize_algorithm_text(text);
    if(normalized.find("bell") != std::string::npos)
    {
        return Bell;
    }
    return AVL;
}

bool parse_bool_text(const std::string& text, bool fallback)
{
    if(text.empty())
    {
        return fallback;
    }
    const std::string normalized = normalize_algorithm_text(text);
    if(normalized == "1" || normalized == "true" || normalized == "yes" || normalized == "on" || normalized == "checked" || normalized == "enable" || normalized == "enabled")
    {
        return true;
    }
    if(normalized == "0" || normalized == "false" || normalized == "no" || normalized == "off" || normalized == "unchecked" || normalized == "disable" || normalized == "disabled")
    {
        return false;
    }
    return fallback;
}

bool parse_strict_bool_text(const std::string& text, bool* value)
{
    if(value == NULL)
    {
        return false;
    }
    const std::string normalized = normalize_algorithm_text(text);
    if(normalized == "1" || normalized == "true" || normalized == "yes" ||
       normalized == "on" || normalized == "checked" || normalized == "enable" ||
       normalized == "enabled")
    {
        *value = true;
        return true;
    }
    if(normalized == "0" || normalized == "false" || normalized == "no" ||
       normalized == "off" || normalized == "unchecked" || normalized == "disable" ||
       normalized == "disabled")
    {
        *value = false;
        return true;
    }
    return false;
}

bool parse_strict_double_text(const std::string& text, double* value)
{
    const char* start = text.c_str();
    char* end = NULL;
    double parsed = 0.0;
    if(value == NULL || start == NULL || *start == '\0')
    {
        return false;
    }
    errno = 0;
    parsed = strtod(start, &end);
    if(end == start || errno == ERANGE || !wmc::is_finite_number(parsed))
    {
        return false;
    }
    while(end != NULL && *end != '\0' && isspace((unsigned char)*end))
    {
        ++end;
    }
    if(end == NULL || *end != '\0')
    {
        return false;
    }
    *value = parsed;
    return true;
}

bool parse_strict_int_text(const std::string& text, int* value)
{
    const char* start = text.c_str();
    char* end = NULL;
    long parsed = 0;
    if(value == NULL || start == NULL || *start == '\0')
    {
        return false;
    }
    errno = 0;
    parsed = strtol(start, &end, 10);
    if(end == start || errno == ERANGE || parsed < INT_MIN || parsed > INT_MAX)
    {
        return false;
    }
    while(end != NULL && *end != '\0' && isspace((unsigned char)*end))
    {
        ++end;
    }
    if(end == NULL || *end != '\0')
    {
        return false;
    }
    *value = (int)parsed;
    return true;
}

double clamp_double_value(double value, double lower, double upper)
{
    if(value < lower)
    {
        return lower;
    }
    if(value > upper)
    {
        return upper;
    }
    return value;
}

double percentile_value(std::vector<double> values, double percentile)
{
    if(values.empty())
    {
        return 0.0;
    }
    if(percentile < 0.0)
    {
        percentile = 0.0;
    }
    if(percentile > 100.0)
    {
        percentile = 100.0;
    }
    std::sort(values.begin(), values.end());
    const double position = (percentile / 100.0) * (double)(values.size() - 1);
    const size_t lower_index = (size_t)std::floor(position);
    const size_t upper_index = (size_t)std::ceil(position);
    if(lower_index == upper_index)
    {
        return values[lower_index];
    }
    const double fraction = position - (double)lower_index;
    return values[lower_index] * (1.0 - fraction) + values[upper_index] * fraction;
}

bool parse_weight_mode_text(const std::string& text, wmc::WeightMode* mode)
{
    const std::string normalized = normalize_algorithm_text(text);
    if(mode == NULL)
    {
        return false;
    }
    if(normalized == "uniform")
    {
        *mode = wmc::WeightMode_Uniform;
        return true;
    }
    if(normalized == "coordinate_gaussian" || normalized == "coordinate gaussian" || normalized == "spatial_gaussian" || normalized == "spatial gaussian" || normalized == "coordinate" || normalized == "coord" || normalized == "coord_gauss")
    {
        *mode = wmc::WeightMode_CoordinateGaussian;
        return true;
    }
    if(normalized == "soft_high_porosity" || normalized == "soft high porosity" ||
       normalized == "soft_high" || normalized == "soft")
    {
        *mode = wmc::WeightMode_SoftHighPorosity;
        return true;
    }
    return false;
}

std::string sanitize_grid_name_base(const std::string& raw)
{
    std::string out;
    bool previous_underscore = false;
    for(size_t i = 0; i < raw.size(); i++)
    {
        const unsigned char c = (unsigned char)raw[i];
        if(isalnum(c))
        {
            out.push_back((char)c);
            previous_underscore = false;
        }
        else if(!previous_underscore)
        {
            out.push_back('_');
            previous_underscore = true;
        }
    }
    while(!out.empty() && out[out.size() - 1] == '_')
    {
        out.erase(out.size() - 1);
    }
    if(out.empty())
    {
        out = "HOSC_MW_OUT";
    }
    const size_t max_name_length = 56;
    if(out.size() > max_name_length)
    {
        unsigned int hash = 2166136261u;
        for(size_t i = 0; i < raw.size(); i++)
        {
            hash ^= (unsigned char)raw[i];
            hash *= 16777619u;
        }
        char hash_buffer[16];
#if defined(_MSC_VER)
        _snprintf(hash_buffer, sizeof(hash_buffer) - 1, "_h%08X", hash);
#else
        snprintf(hash_buffer, sizeof(hash_buffer) - 1, "_h%08X", hash);
#endif
        hash_buffer[sizeof(hash_buffer) - 1] = '\0';
        const std::string hash_suffix(hash_buffer);
        const size_t prefix_length = max_name_length > hash_suffix.size()
            ? max_name_length - hash_suffix.size()
            : max_name_length;
        out.erase(prefix_length);
        while(!out.empty() && out[out.size() - 1] == '_')
        {
            out.erase(out.size() - 1);
        }
        out += hash_suffix;
    }
    return out;
}

std::string compact_double_token(double value)
{
    char buffer[64];
#if defined(_MSC_VER)
    _snprintf(buffer, sizeof(buffer) - 1, "%.4g", value);
#else
    snprintf(buffer, sizeof(buffer) - 1, "%.4g", value);
#endif
    buffer[sizeof(buffer) - 1] = '\0';

    std::string out;
    for(size_t i = 0; buffer[i] != '\0'; i++)
    {
        const char c = buffer[i];
        if(c == '-')
        {
            out += "m";
        }
        else if(c == '+')
        {
            continue;
        }
        else if(c == '.')
        {
            out += "p";
        }
        else if(isalnum((unsigned char)c))
        {
            out.push_back(c);
        }
    }
    return out.empty() ? "0" : out;
}

std::string shape_token(const GridShape* shape)
{
    if(shape == NULL)
    {
        return "0x0x0";
    }
    std::ostringstream oss;
    oss << shape->nx << "x" << shape->ny << "x" << shape->nz;
    return oss.str();
}

const char* target_name_token(int target)
{
    switch(target)
    {
    case 3: return "MOM";
    case 4: return "CUM";
    default: return "OUT";
    }
}

const char* algorithm_name_token(Algorithm algorithm)
{
    return algorithm == Bell ? "BELL" : "AVL";
}

const char* weight_name_token(wmc::WeightMode mode)
{
    switch(mode)
    {
    case wmc::WeightMode_CoordinateGaussian: return "WCG";
    case wmc::WeightMode_SoftHighPorosity: return "WSH";
    case wmc::WeightMode_Uniform:
    default:
        return "WU";
    }
}

bool safe_get_source_data(HOSC_WINDOW* self, const Parameters_handler* parameters, Error_messages_handler* error_mesgs, unsigned long* exception_code)
{
    if(exception_code != NULL)
    {
        *exception_code = 0;
    }
    __try
    {
        return self->get_Source_Data(parameters, error_mesgs);
    }
    __except(EXCEPTION_EXECUTE_HANDLER)
    {
        if(exception_code != NULL)
        {
            *exception_code = 1;
        }
        return false;
    }
}

bool safe_get_calculate_algorithm(HOSC_WINDOW* self, const Parameters_handler* parameters, Error_messages_handler* error_mesgs, unsigned long* exception_code)
{
    if(exception_code != NULL)
    {
        *exception_code = 0;
    }
    __try
    {
        append_hosc_window_trace("safe_get_Calculate_Algorithm.before_call");
        const bool ok = self->get_Calculate_Algorithm(parameters, error_mesgs);
        append_hosc_window_tracef("safe_get_Calculate_Algorithm.after_call ok=%d algorithm=%d", ok ? 1 : 0, (int)self->al);
        return ok;
    }
    __except(EXCEPTION_EXECUTE_HANDLER)
    {
        if(exception_code != NULL)
        {
            *exception_code = 1;
        }
        self->al = AVL;
        append_hosc_window_trace("safe_get_Calculate_Algorithm.exception fallback=AVL");
        return true;
    }
}

bool safe_get_weighted_moment_options(HOSC_WINDOW* self, const Parameters_handler* parameters, Error_messages_handler* error_mesgs, unsigned long* exception_code)
{
    if(exception_code != NULL)
    {
        *exception_code = 0;
    }
    __try
    {
        append_hosc_window_trace("safe_get_Weighted_Moment_Options.before_call");
        const bool ok = self->get_Weighted_Moment_Options(parameters, error_mesgs);
        append_hosc_window_tracef(
            "safe_get_Weighted_Moment_Options.after_call ok=%d enabled=%d mode=%d",
            ok ? 1 : 0,
            self->useWeightedMoment ? 1 : 0,
            (int)self->weightedOptions.mode);
        return ok;
    }
    __except(EXCEPTION_EXECUTE_HANDLER)
    {
        if(exception_code != NULL)
        {
            *exception_code = 1;
        }
        self->useWeightedMoment = false;
        self->weightedOptions = wmc::WeightedMomentOptions();
        append_hosc_window_trace("safe_get_Weighted_Moment_Options.exception fallback=disabled");
        return true;
    }
}

bool safe_get_windows_moment(HOSC_WINDOW* self, const Parameters_handler* parameters, Error_messages_handler* error_mesgs, unsigned long* exception_code)
{
    if(exception_code != NULL)
    {
        *exception_code = 0;
    }
    __try
    {
        append_hosc_window_trace("safe_get_Windows_Moment.before_call");
        const bool ok = self->get_Windows_Moment(parameters, error_mesgs);
        append_hosc_window_tracef("safe_get_Windows_Moment.after_call ok=%d", ok ? 1 : 0);
        return ok;
    }
    __except(EXCEPTION_EXECUTE_HANDLER)
    {
        const unsigned long code = (unsigned long)GetExceptionCode;
        if(exception_code != NULL)
        {
            *exception_code = code;
        }
        append_hosc_window_tracef("safe_get_Windows_Moment.exception code=0x%08lx", code);
        return false;
    }
}

int safe_calculate_moment_window(HOSC_WINDOW* self, unsigned long* exception_code)
{
    if(exception_code != NULL)
    {
        *exception_code = 0;
    }
    __try
    {
        return self->calculateMomentWindow();
    }
    // A structured exception can follow heap corruption.  Do not turn it into
    // a normal plugin failure and leave the host process running.  The audit
    // directory remains in its last durable state (normally "started") for
    // the debugger/PDB-assisted investigation.
    __except(EXCEPTION_CONTINUE_SEARCH)
    {
        return -1;
    }
}

bool safe_project_new_object(GsTL_project* project, const std::string* object_name, unsigned long* exception_code)
{
    if(exception_code != NULL)
    {
        *exception_code = 0;
    }
    if(project == NULL || object_name == NULL || object_name->empty())
    {
        return false;
    }
    try
    {
        project->new_object(*object_name);
        return true;
    }
    catch(const std::exception& e)
    {
        append_hosc_window_tracef("safe_project_new_object.exception message=%s", e.what());
        if(exception_code != NULL)
        {
            *exception_code = 1;
        }
        return false;
    }
    catch(...)
    {
        append_hosc_window_trace("safe_project_new_object.unknown_exception");
        if(exception_code != NULL)
        {
            *exception_code = 1;
        }
        return false;
    }
}

bool project_execute_action_unchecked(
    GsTL_project* project,
    const std::string& action_name,
    const std::string& params,
    Error_messages_handler* errors)
{
    return project->execute(action_name, params, errors);
}

bool safe_project_execute_action(
    GsTL_project* project,
    const std::string* action_name,
    const std::string* params,
    Error_messages_handler* errors,
    unsigned long* exception_code)
{
    if(exception_code != NULL)
    {
        *exception_code = 0;
    }
    if(project == NULL || action_name == NULL || action_name->empty() || params == NULL)
    {
        return false;
    }
    __try
    {
        return project_execute_action_unchecked(project, *action_name, *params, errors);
    }
    __except(EXCEPTION_EXECUTE_HANDLER)
    {
        if(exception_code != NULL)
        {
            *exception_code = 1;
        }
        append_hosc_window_trace("safe_project_execute_action.exception");
        return false;
    }
}

double bell_binomial(int n, int k)
{
    if(k < 0 || k > n)
    {
        return 0.0;
    }
    if(k == 0 || k == n)
    {
        return 1.0;
    }
    if(k > n - k)
    {
        k = n - k;
    }

    long double result = 1.0;
    for(int i = 1; i <= k; i++)
    {
        result *= (long double)(n - k + i);
        result /= (long double)i;
    }
    return (double)result;
}

double calculate_window_bell_cumulant_from_samples(const std::vector<double>& samples, int order)
{
    if(order <= 0 || samples.empty())
    {
        return 0.0;
    }

    std::vector<double> moments(order + 1, 0.0);
    std::vector<double> cumulants(order + 1, 0.0);
    for(size_t i = 0; i < samples.size(); i++)
    {
        long double power = 1.0;
        for(int n = 1; n <= order; n++)
        {
            power *= (long double)samples[i];
            moments[n] += (double)power;
        }
    }

    for(int n = 1; n <= order; n++)
    {
        moments[n] /= (double)samples.size();
    }

    cumulants[1] = moments[1];
    for(int n = 2; n <= order; n++)
    {
        long double value = moments[n];
        for(int j = 1; j < n; j++)
        {
            value -= (long double)bell_binomial(n - 1, j - 1) * (long double)cumulants[j] * (long double)moments[n - j];
        }
        cumulants[n] = (double)value;
    }

    return cumulants[order];
}

double calculate_window_moment_from_samples(const std::vector<double>& samples, int order)
{
    if(samples.empty())
    {
        return 0.0;
    }
    if(order == 0)
    {
        return 1.0;
    }
    long double sum = 0.0;
    for(size_t i = 0; i < samples.size(); i++)
    {
        sum += pow(samples[i], (double)order);
    }
    return (double)(sum / (long double)samples.size());
}

bool calculate_window_avl_cumulant_from_samples(const std::vector<double>& samples, int order, double* result);

bool calculate_window_cumulant_from_samples(
	const std::vector<double>& samples, int order, Algorithm algorithm, double* result)
{
	if(result == NULL)
	{
		return false;
	}
	if(algorithm == Bell)
	{
		*result = calculate_window_bell_cumulant_from_samples(samples, order);
		return wmc::is_finite_number(*result);
	}
	return calculate_window_avl_cumulant_from_samples(samples, order, result);
}

std::vector<double> centered_unweighted_samples(const std::vector<double>& samples)
{
	std::vector<double> centered(samples);
	long double sum = 0.0;
	if(samples.empty())
	{
		return centered;
	}
	for(size_t i = 0; i < samples.size(); i++)
	{
		sum += (long double)samples[i];
	}
	const double mean = (double)(sum / (long double)samples.size());
	for(size_t i = 0; i < centered.size(); i++)
	{
		centered[i] -= mean;
	}
	return centered;
}

bool is_float_representable(double value)
{
	const double limit = (double)std::numeric_limits<float>::max();
	return wmc::is_finite_number(value) && value >= -limit && value <= limit;
}

const char* const HOSC_WINDOW_BUILD_ID = "hosc_window_20261020_r7_repeat_diagnostic";

bool create_window_audit_directory(
	const std::string& base_directory,
	unsigned long process_id,
	unsigned long run_sequence,
	std::string* run_directory)
{
	if(run_directory == NULL || base_directory.empty())
	{
		return false;
	}
	if(_mkdir(base_directory.c_str()) != 0 && errno != EEXIST)
	{
		return false;
	}
	char suffix[160];
	const long now = (long)time(NULL);
	const unsigned long tick = (unsigned long)clock();
	_snprintf(
		suffix,
		sizeof(suffix) - 1,
		"HOSC_WINDOW_%ld_pid%lu_seq%lu_tick%lu",
		now,
		process_id,
		run_sequence,
		tick);
	suffix[sizeof(suffix) - 1] = '\0';
	*run_directory = base_directory;
	if(!run_directory->empty() && (*run_directory)[run_directory->size() - 1] != '\\' && (*run_directory)[run_directory->size() - 1] != '/')
	{
		*run_directory += "\\";
	}
	*run_directory += suffix;
	if(_mkdir(run_directory->c_str()) != 0)
	{
		return false;
	}
	return true;
}

bool write_window_audit_status(
	const std::string& run_directory,
	const std::string& run_id,
	unsigned long process_id,
	unsigned long run_sequence,
	const char* state,
	const char* detail)
{
	if(run_directory.empty() || state == NULL)
	{
		return false;
	}
	std::ofstream status((run_directory + "\\run_status.txt").c_str(), std::ios::out | std::ios::trunc);
	if(!status)
	{
		return false;
	}
	status << "schema_version=hosc_window_run_status.v1\n";
	status << "build_id=" << HOSC_WINDOW_BUILD_ID << "\n";
	status << "run_id=" << run_id << "\n";
	status << "process_id=" << process_id << "\n";
	status << "run_sequence=" << run_sequence << "\n";
	status << "state=" << state << "\n";
	status << "detail=" << (detail == NULL ? "" : detail) << "\n";
	status << "updated_epoch=" << (long)time(NULL) << "\n";
	return status.good();
}

bool write_window_audit_sidecar(
	const std::string& run_directory,
	const std::string& output_name,
	const std::string& source_grid_name,
	const std::string& source_property_name,
	int target,
	Algorithm algorithm,
	bool weighted,
	const wmc::WeightedMomentOptions& weights,
	bool center_cumulant,
	bool output_standardized,
	bool output_diagnostics,
	bool invert_output,
	bool enable_percentile_clipping,
	double percentile_vmin,
	double percentile_vmax,
	double variance_epsilon,
	const std::string& run_id,
	unsigned long process_id,
	unsigned long run_sequence,
	int order,
	int window_x,
	int window_y,
	int window_z,
	int step_x,
	int step_y,
	int step_z,
	double origin_x,
	double origin_y,
	double origin_z,
	double spacing_x,
	double spacing_y,
	double spacing_z,
	int out_x,
	int out_y,
	int out_z,
	const std::vector<double>& values,
	const std::vector<double>& support,
	const std::vector<double>& ess,
	const std::vector<double>& numeric_valid,
	const std::vector<double>& display_valid,
	const std::vector<double>& kappa2,
	const std::vector<double>& lambda,
	const std::vector<double>& lambda_valid)
{
	if(run_directory.empty())
	{
		return false;
	}
	std::ofstream config((run_directory + "\\resolved_config.txt").c_str(), std::ios::out | std::ios::trunc);
	if(!config)
	{
		return false;
	}
	config.precision(17);
	config << "schema_version=hosc_window_audit.v3\n";
	config << "build_id=" << HOSC_WINDOW_BUILD_ID << "\n";
	config << "run_id=" << run_id << "\n";
	config << "process_id=" << process_id << "\n";
	config << "run_sequence=" << run_sequence << "\n";
	config << "weighting_formula_version=" << wmc::weighting_formula_version() << "\n";
	config << "source_grid=" << source_grid_name << "\n";
	config << "source_property=" << source_property_name << "\n";
	config << "output_property=" << output_name << "\n";
	config << "target=" << (target == 4 ? "Cumulant" : "Moment") << "\n";
	config << "unweighted_algorithm=" << (algorithm == Bell ? "Bell" : "AVL") << "\n";
	config << "weighted_route=" << (weighted ? "shared_weighted_recurrence" : "none") << "\n";
	config << "weight_mode=" << wmc::weight_mode_name(weights.mode) << "\n";
	config << "normalize_weights=" << (weights.normalize_weights ? "true" : "false") << "\n";
	config << "weight_bandwidth=" << weights.bandwidth << "\n";
	config << "coordinate_center_mode=" << (weights.mode == wmc::WeightMode_CoordinateGaussian ? "window_center" : "not_applicable") << "\n";
	config << "high_porosity_threshold=" << weights.high_porosity_threshold << "\n";
	config << "soft_transition_width=" << weights.soft_transition_width << "\n";
	config << "soft_weight_strength=" << weights.soft_weight_strength << "\n";
	config << "min_weight=" << weights.min_weight << "\n";
	config << "zero_first_order_moment=" << (weights.zero_first_order_moment ? "true" : "false") << "\n";
	config << "center_cumulant_input_for_order_ge_2=" << (center_cumulant ? "true" : "false") << "\n";
	config << "output_standardized_cumulant=" << (output_standardized ? "true" : "false") << "\n";
	config << "output_diagnostics=" << (output_diagnostics ? "true" : "false") << "\n";
	config << "invert_window_output=" << (invert_output ? "true" : "false") << "\n";
	config << "enable_percentile_clipping=" << (enable_percentile_clipping ? "true" : "false") << "\n";
	config << "percentile_vmin=" << percentile_vmin << "\n";
	config << "percentile_vmax=" << percentile_vmax << "\n";
	config << "variance_epsilon=" << variance_epsilon << "\n";
	config << "moment_order=" << order << "\n";
	config << "window=" << window_x << "," << window_y << "," << window_z << "\n";
	config << "step=" << step_x << "," << step_y << "," << step_z << "\n";
	config << "source_origin=" << origin_x << "," << origin_y << "," << origin_z << "\n";
	config << "source_spacing=" << spacing_x << "," << spacing_y << "," << spacing_z << "\n";
	config << "output_dimensions=" << out_x << "," << out_y << "," << out_z << "\n";
	config.close();

	std::ofstream rows((run_directory + "\\window_results.csv").c_str(), std::ios::out | std::ios::trunc);
	if(!rows)
	{
		return false;
	}
	rows.precision(17);
	rows << "node_id,start_x,start_y,start_z,center_x,center_y,center_z,raw_value,kappa2,lambda,support,ess,numeric_valid,lambda_requested,lambda_valid,display_valid,reason\n";
	const size_t row_count = values.size();
	for(size_t index = 0; index < row_count; ++index)
	{
		const int node = (int)index;
		const int ox = node % out_x;
		const int oy = (node / out_x) % out_y;
		const int oz = node / (out_x * out_y);
		const int start_x = ox * step_x;
		const int start_y = oy * step_y;
		const int start_z = oz * step_z;
		const double center_x = origin_x + ((double)start_x + 0.5 * (double)(window_x - 1)) * spacing_x;
		const double center_y = origin_y + ((double)start_y + 0.5 * (double)(window_y - 1)) * spacing_y;
		const double center_z = origin_z + ((double)start_z + 0.5 * (double)(window_z - 1)) * spacing_z;
		const bool numeric = index < numeric_valid.size() && numeric_valid[index] != 0.0;
		const bool displayed = index < display_valid.size() && display_valid[index] != 0.0;
		const bool standardized = index < lambda_valid.size() && lambda_valid[index] != 0.0;
		const bool lambda_requested = output_standardized && target == 4;
		const char* reason = !numeric ? "incomplete_or_invalid_support" :
			(!displayed ? "float_overflow" :
			 (!lambda_requested ? "lambda_not_requested" :
			  (!standardized ? "variance_below_epsilon_or_lambda_nonfinite" : "ok")));
		rows << node << "," << start_x << "," << start_y << "," << start_z << ","
			 << center_x << "," << center_y << "," << center_z << ","
			 << values[index] << ","
			 << (index < kappa2.size() ? kappa2[index] : 0.0) << ","
			 << (index < lambda.size() ? lambda[index] : 0.0) << ","
			 << (index < support.size() ? support[index] : 0.0) << ","
			 << (index < ess.size() ? ess[index] : 0.0) << ","
			 << (numeric ? 1 : 0) << "," << (lambda_requested ? 1 : 0) << "," << (standardized ? 1 : 0) << "," << (displayed ? 1 : 0) << "," << reason << "\n";
	}
	return rows.good();
}

void reset_window_moment_workspace()
{
    if(G_Moment != NULL)
    {
        free(G_Moment);
        G_Moment = NULL;
    }
    if(G_MomentIndex != NULL)
    {
        if(G_SeriesNum > 0)
        {
            Free2DVec_char(&G_MomentIndex, G_SeriesNum);
        }
        else
        {
            free(G_MomentIndex);
            G_MomentIndex = NULL;
        }
    }
    if(G_MomentHash != NULL)
    {
        free(G_MomentHash);
        G_MomentHash = NULL;
    }
    G_MomentHashSize = 0;
    G_MomentHashReady = false;
    G_MomentFull = false;
    G_SeriesNum = 0;
    G_CurrMomentPos = 0;
}

void clear_window_avl_roots(int order)
{
    if(RootAVL2 == NULL)
    {
        return;
    }
    for(int i = 1; i <= order; i++)
    {
        if(RootAVL2[i] != NULL)
        {
            FreeAVLTree(&RootAVL2[i]);
            RootAVL2[i] = NULL;
        }
    }
    free(RootAVL2);
    RootAVL2 = NULL;
}

void free_window_template_samples(int sample_count)
{
    if(G_Template == NULL)
    {
        return;
    }
    for(int i = 0; i < sample_count; i++)
    {
        if(G_Template[i].Value != NULL)
        {
            free(G_Template[i].Value);
            G_Template[i].Value = NULL;
        }
    }
    if(G_Template->index != NULL)
    {
        free(G_Template->index);
        G_Template->index = NULL;
    }
    free(G_Template);
    G_Template = NULL;
    G_TemplateNum = 0;
}

bool calculate_window_avl_cumulant_from_samples(const std::vector<double>& samples, int order, double* result)
{
    if(result == NULL || samples.empty() || order <= 0)
    {
        return false;
    }

    G_Template = (Template*)calloc(samples.size(), sizeof(Template));
    if(G_Template == NULL)
    {
        return false;
    }
    G_TemplateNum = (int)samples.size();

    for(size_t i = 0; i < samples.size(); i++)
    {
        G_Template[i].Value = (double*)calloc(2, sizeof(double));
        if(G_Template[i].Value == NULL)
        {
            free_window_template_samples((int)samples.size());
            return false;
        }
        G_Template[i].Value[1] = samples[i];
    }

    G_Template->index = (int*)calloc(2, sizeof(int));
    if(G_Template->index == NULL)
    {
        free_window_template_samples((int)samples.size());
        return false;
    }
    G_Template->index[0] = order;
    G_Template->index[1] = order;

    reset_window_moment_workspace();
    clear_window_avl_roots(1);
    if(InitRoot2(1) == -1)
    {
        free_window_template_samples((int)samples.size());
        return false;
    }

    *result = BetterRecCumulant3(G_Template->index, 1);
    RefreshTree(1);
    clear_window_avl_roots(1);
    reset_window_moment_workspace();
    free_window_template_samples((int)samples.size());
    return true;
}

bool is_effective_2d_grid(const RGrid* grid_obj)
{
    if(grid_obj == NULL)
    {
        return true;
    }
	return grid_obj->nz() <= 1;
}

bool is_effective_2d_dims(int nx, int ny, int nz)
{
	(void)nx;
	(void)ny;
	return nz <= 1;
}

int round_grid_coordinate(double value)
{
    if(value >= 0.0)
        return (int)floor(value + 0.5);
    return (int)ceil(value - 0.5);
}

void write_be_uint32(std::ofstream& out, unsigned int value)
{
    char bytes[4];
    bytes[0] = (char)((value >> 24) & 0xFF);
    bytes[1] = (char)((value >> 16) & 0xFF);
    bytes[2] = (char)((value >> 8) & 0xFF);
    bytes[3] = (char)(value & 0xFF);
    out.write(bytes, 4);
}

void write_be_float(std::ofstream& out, float value)
{
    union
    {
        float f;
        unsigned int u;
    } bits;
    bits.f = value;
    write_be_uint32(out, bits.u);
}

void write_qdatastream_string(std::ofstream& out, const std::string& value)
{
    write_be_uint32(out, (unsigned int)value.size() + 1U);
    out.write(value.c_str(), (std::streamsize)value.size() + 1);
}

bool write_sgems_pointset_file(
    const std::string& filename,
    const std::string& object_name,
    const std::string& property_name,
    const std::vector<Geostat_grid::location_type>& point_locations,
    const std::vector<double>& point_values)
{
    if(point_locations.size() != point_values.size())
    {
        append_hosc_window_tracef(
            "write_sgems_pointset_file.fail size_mismatch locations=%d values=%d",
            (int)point_locations.size(),
            (int)point_values.size());
        return false;
    }

    std::ofstream out(filename.c_str(), std::ios::binary);
    if(!out)
    {
        append_hosc_window_tracef("write_sgems_pointset_file.fail open path=%s", filename.c_str());
        return false;
    }

    write_be_uint32(out, 0xB211175D);
    write_qdatastream_string(out, "Point_set");
    write_qdatastream_string(out, object_name);
    write_be_uint32(out, 100);
    write_be_uint32(out, (unsigned int)point_locations.size());
    write_be_uint32(out, 1);
    write_qdatastream_string(out, property_name);

    for(size_t i = 0; i < point_locations.size(); i++)
    {
        write_be_float(out, (float)point_locations[i].x());
        write_be_float(out, (float)point_locations[i].y());
        write_be_float(out, (float)point_locations[i].z());
    }
    for(size_t i = 0; i < point_values.size(); i++)
    {
        write_be_float(out, (float)point_values[i]);
    }

    if(!out)
    {
        append_hosc_window_tracef("write_sgems_pointset_file.fail write path=%s", filename.c_str());
        return false;
    }
    append_hosc_window_tracef(
        "write_sgems_pointset_file.end path=%s name=%s points=%d property=%s",
        filename.c_str(),
        object_name.c_str(),
        (int)point_locations.size(),
        property_name.c_str());
    return true;
}
}










HOSC_WINDOW::HOSC_WINDOW()
{
    source_grid_ = 0;
	source_rgrid_ = 0;
    source_x = 0;
	source_y = 0;
	source_z = 0;
	source_val = 0;
    target_grid_ = 0;

	grid = NULL;
	grid_alloc_rows_ = 0;
	source_nx_ = 0;
	source_ny_ = 0;
	source_nz_ = 0;
	source_grid_name_.clear();
	source_property_name_.clear();
	source_min_x_ = 1;
	source_min_y_ = 1;
	source_min_z_ = 1;
	source_present_.clear();
	Key = NULL;
	al = AVL;
	current_project_ = NULL;
	auto_grid_ = NULL;
	auto_property_ = NULL;
	auto_grid_name_.clear();
	auto_grid_registered_ = false;
	output_keep_mask_.clear();
	target_grid_ = NULL;
	target_property_ = NULL;
	momentWindowL = 0;
	momentWindowW = 0;
	momentWindowH = 0;
	momentOrder = 0;
	stepX = 0;
	stepY = 0;
	stepZ = 0;
	invertWindowOutput = false;
	enablePercentileClipping = false;
	outputDiagnostics = true;
	centerCumulantInputForOrderGE2 = true;
	outputStandardizedCumulant = false;
	writeAuditSidecar = false;
	percentileClipMin = 2.0;
	percentileClipMax = 98.0;
	varianceEpsilon = 1.0e-12;
	auditOutputDirectory.clear();
	auditRunDirectory.clear();
	runId_.clear();
	processId_ = (unsigned long)GetCurrentProcessId();
	runSequence_ = 0;
	useWeightedMoment = false;
	coordinateWeightUsesWindowCenter = true;
	weightedOptions = wmc::WeightedMomentOptions();
	// Fixed public policy; legacy PAR fields below are checked, not applied.
	weightedOptions.normalize_weights = true;
	weightedOptions.zero_first_order_moment = false;

	G_CalculateTarget = -1;
	G_TargetDimension = 0;
	append_hosc_window_tracef(
		"HOSC_WINDOW.ctor build=%s instance=%p",
		HOSC_WINDOW_BUILD_ID,
		this);
}

HOSC_WINDOW::~HOSC_WINDOW()
{
	append_hosc_window_instance_trace(this, runSequence_, "dtor.body.begin");
	append_hosc_window_tracef(
		"HOSC_WINDOW.dtor.resources instance=%p run_sequence=%lu grid=%p rows=%d source_grid=%p auto_grid=%p auto_property=%p",
		this,
		runSequence_,
		grid,
		grid_alloc_rows_,
		source_grid_,
		auto_grid_,
		auto_property_);
	release_source_grid_cache();
	source_grid_ = NULL;
	source_rgrid_ = NULL;
	source_x = NULL;
	source_y = NULL;
	source_z = NULL;
	source_val = NULL;
	target_grid_ = NULL;
	target_property_ = NULL;
	auto_grid_ = NULL;
	auto_property_ = NULL;
	auto_grid_name_.clear();
	auto_grid_registered_ = false;
	current_project_ = NULL;
	Key = NULL;
	append_hosc_window_instance_trace(this, runSequence_, "dtor.body.end.members_and_base_destruction_next");
}

std::string HOSC_WINDOW::build_output_grid_base_name(const GridShape* tempGridShape, bool pointset) const
{
	const char target_token = G_CalculateTarget == 4 ? 'C' : (G_CalculateTarget == 3 ? 'M' : 'O');
	const char algorithm_token = al == Bell ? 'B' : 'A';
	char weight_token = 'U';
	if(useWeightedMoment)
	{
		if(weightedOptions.mode == wmc::WeightMode_CoordinateGaussian)
		{
			weight_token = 'C';
		}
		else if(weightedOptions.mode == wmc::WeightMode_SoftHighPorosity)
		{
			weight_token = 'S';
		}
		else
		{
			weight_token = 'W';
		}
	}
	const int output_dim = (tempGridShape != NULL && tempGridShape->nx > 1 && tempGridShape->ny > 1 && tempGridShape->nz > 1) ? 3 : 2;

	std::ostringstream name;
	name << "myWM"
	     << target_token
	     << momentOrder
	     << algorithm_token
	     << weight_token
	     << output_dim;
	if(pointset)
	{
		name << "P";
	}
	if(invertWindowOutput)
	{
		name << "I";
	}
	if(enablePercentileClipping)
	{
		name << "Q";
	}
	return sanitize_grid_name_base(name.str());
}

void HOSC_WINDOW::release_source_grid_cache()
{
	if(grid == NULL)
	{
		grid_alloc_rows_ = 0;
		return;
	}
	append_hosc_window_tracef("release_source_grid_cache.begin grid=%p rows=%d", grid, grid_alloc_rows_);
	for(int i = 0; i < grid_alloc_rows_; i++)
	{
		if(grid[i] != NULL)
		{
			delete [] grid[i];
			grid[i] = NULL;
		}
	}
	delete [] grid;
	grid = NULL;
	grid_alloc_rows_ = 0;
	append_hosc_window_trace("release_source_grid_cache.end");
}


bool HOSC_WINDOW::initialize( const Parameters_handler* parameters,
							   Error_messages_handler* error_mesgs )
{
	hosc_validation::FeedbackGuard feedback(error_mesgs, "hosc_window_validation_summary");
	processId_ = (unsigned long)GetCurrentProcessId();
	runSequence_ = next_hosc_window_run_sequence();
	char run_id_buffer[96];
	_snprintf(
		run_id_buffer,
		sizeof(run_id_buffer) - 1,
		"pid%lu_seq%lu",
		processId_,
		runSequence_);
	run_id_buffer[sizeof(run_id_buffer) - 1] = '\0';
	runId_ = run_id_buffer;
	append_hosc_window_instance_trace(this, runSequence_, "initialize.begin");
	GsTLlog << "HOSC_WINDOW initialize start" << gstlIO::end;
	if(!get_Calculate_Target(parameters, error_mesgs))
	{
		append_hosc_window_trace("initialize.fail target");
		return false;
	}
	append_hosc_window_tracef("initialize.after_target target=%d", G_CalculateTarget);
	if(!get_Calculate_Algorithm(parameters, error_mesgs))
	{
		append_hosc_window_trace("initialize.fail algorithm");
		return false;
	}
	append_hosc_window_trace("initialize.after_algorithm");
	if(!get_Weighted_Moment_Options(parameters, error_mesgs))
	{
		append_hosc_window_trace("initialize.fail weighted_options");
		return false;
	}
	append_hosc_window_tracef(
		"initialize.after_weighted enabled=%d mode=%d bandwidth=%.17g reference=%.17g min=%.17g normalize=%d zero_first=%d",
		useWeightedMoment ? 1 : 0,
		(int)weightedOptions.mode,
		weightedOptions.bandwidth,
		weightedOptions.reference_value,
		weightedOptions.min_weight,
		weightedOptions.normalize_weights ? 1 : 0,
		weightedOptions.zero_first_order_moment ? 1 : 0);

	append_hosc_window_trace("initialize.before_source_data");
    if ( !get_Source_Data(parameters, error_mesgs) )
	{
		append_hosc_window_trace("initialize.fail source_data");
        return false;
	}
	append_hosc_window_trace("initialize.after_source_data");
	append_hosc_window_instance_trace(this, runSequence_, "source_ready.host_log.begin");
	GsTLlog << "HOSC_WINDOW source ready: nx=" << source_nx_ << ", ny=" << source_ny_ << ", nz=" << source_nz_ << gstlIO::end;
	append_hosc_window_instance_trace(this, runSequence_, "source_ready.host_log.end");
	append_hosc_window_tracef(
		"initialize.source_ready nx=%d ny=%d nz=%d nxyz=%d raw_nxyz=%d grid=%p source_val=%p",
		source_nx_,
		source_ny_,
		source_nz_,
		(int)grid[0][0],
		source_grid_->size(),
		source_grid_,
		source_val);











	append_hosc_window_instance_trace(this, runSequence_, "window_params.call.begin");
	unsigned long window_parameters_exception = 0;
	const bool window_parameters_ok = safe_get_windows_moment(
		this,
		parameters,
		error_mesgs,
		&window_parameters_exception);
	append_hosc_window_tracef(
		"run_sequence=%lu instance=%p stage=window_params.call.returned ok=%d exception=0x%08lx",
		runSequence_,
		this,
		window_parameters_ok ? 1 : 0,
		window_parameters_exception);
	if(!window_parameters_ok)
	{
		append_hosc_window_trace("initialize.fail windows_moment");
		return false;
	}
	if(error_mesgs == NULL)
	{
		append_hosc_window_trace("initialize.fail error_messages_null");
		return false;
	}
	append_hosc_window_tracef(
		"initialize.after_window_params error_messages_empty=%d",
		error_mesgs->empty() ? 1 : 0);
	append_hosc_window_instance_trace(this, runSequence_, "window_ready.host_log.begin");
	GsTLlog << "HOSC_WINDOW window ready: wx=" << momentWindowL << ", wy=" << momentWindowW << ", wz=" << momentWindowH << ", order=" << momentOrder << ", out_nx=" << myWindowShape.nx << ", out_ny=" << myWindowShape.ny << ", out_nz=" << myWindowShape.nz << gstlIO::end;
	append_hosc_window_instance_trace(this, runSequence_, "window_ready.host_log.end");
	append_hosc_window_tracef(
		"initialize.window_ready order=%d window=%d,%d,%d step=%d,%d,%d out=%d,%d,%d",
		momentOrder,
		momentWindowL,
		momentWindowW,
		momentWindowH,
		stepX,
		stepY,
		stepZ,
		myWindowShape.nx,
		myWindowShape.ny,
		myWindowShape.nz);















	if( !error_mesgs->empty() )
	{
		append_hosc_window_trace("initialize.fail error_messages_not_empty");
		return false;
	}


	GsTLlog << "HOSC_WINDOW initialize success" << gstlIO::end;
	append_hosc_window_trace("initialize.end success");

	return true;
}
int HOSC_WINDOW::execute( GsTL_project* proj )
{
	GsTLlog << "HOSC_WINDOW execute start" << gstlIO::end;
	append_hosc_window_tracef("execute.begin proj=%p source_grid=%p grid=%p auto_grid=%p", proj, source_grid_, grid, auto_grid_);

    if(source_rgrid_ != NULL)
    {
        source_rgrid_->set_level(1);
	    GsTLlog << "HOSC_WINDOW execute source grid level set" << gstlIO::end;
	    append_hosc_window_trace("execute.after_set_level");
    }
    else
    {
        append_hosc_window_trace("execute.skip_set_level_non_regular_source");
    }



	if(grid == NULL)
	{

		append_hosc_window_trace("execute.fail grid_null");
		if(writeAuditSidecar && !auditRunDirectory.empty())
		{
			write_window_audit_status(
				auditRunDirectory,
				runId_,
				processId_,
				runSequence_,
				"failed",
				"execute_grid_null");
		}
		return 1;
	}




#if 0
	int totalDim = Key->NDIR;
	for(int i = 0; i < Key->NDIR; i++)
	{
		if(Key->NLAG[i] <= 1)
			totalDim--;
	}
#endif
	if(G_CalculateTarget == 3 || G_CalculateTarget == 4)
	{
		current_project_ = proj;
		auto_grid_ = NULL;
		auto_property_ = NULL;
		auto_grid_name_.clear();
		auto_grid_registered_ = false;
		output_keep_mask_.clear();
		append_hosc_window_trace("execute.sparse_output_deferred_to_calculate");
		unsigned long calc_exception = 0;
		int calc_status = safe_calculate_moment_window(this, &calc_exception);
		if(calc_status == -1)
		{
			GsTLcerr << "failed to calculate moment window" << gstlIO::end;
			append_hosc_window_tracef("execute.fail calculateMomentWindow exception=0x%08lx", calc_exception);
			if(writeAuditSidecar && !auditRunDirectory.empty())
			{
				write_window_audit_status(
					auditRunDirectory,
					runId_,
					processId_,
					runSequence_,
					"failed",
					"calculate_moment_window_failed");
			}
			return -1;
		}
		if(current_project_ != NULL && !auto_grid_name_.empty() && !auto_grid_registered_)
		{
			append_hosc_window_tracef("execute.before_project_new_object name=%s", auto_grid_name_.c_str());
			unsigned long project_exception = 0;
			if(safe_project_new_object(current_project_, &auto_grid_name_, &project_exception))
			{
				append_hosc_window_trace("execute.after_project_new_object");
			}
			else
			{
				append_hosc_window_tracef("execute.warn project_new_object_failed exception=0x%08lx", project_exception);
			}
		}
		append_hosc_window_trace("execute.after_calculateMomentWindow");
	}

	append_hosc_window_trace("execute.end success");
    return 0;
}

bool HOSC_WINDOW::get_Source_Data(const Parameters_handler* parameters,
                                   Error_messages_handler* error_mesgs )
{
    append_hosc_window_trace("get_Source_Data.enter sparse-fill");
    release_source_grid_cache();
    source_grid_ = NULL;
    source_rgrid_ = NULL;
    source_x = NULL;
    source_y = NULL;
    source_z = NULL;
    source_val = NULL;
    source_nx_ = 0;
    source_ny_ = 0;
    source_nz_ = 0;
	source_grid_name_.clear();
	source_property_name_.clear();
    source_min_x_ = 1;
    source_min_y_ = 1;
    source_min_z_ = 1;
    source_origin_x_ = 1.0;
    source_origin_y_ = 1.0;
    source_origin_z_ = 1.0;
    source_spacing_x_ = 1.0;
    source_spacing_y_ = 1.0;
    source_spacing_z_ = 1.0;
    source_present_.clear();

    string grid_name = parameter_value_compat(parameters, "Source_Grid.value");
    append_hosc_window_tracef("get_Source_Data.begin grid_name=%s", grid_name.c_str());
    if(grid_name.empty())
    {
        error_mesgs->report("Source_Grid", "No source grid selected");
        append_hosc_window_trace("get_Source_Data.fail empty_grid_name");
        return false;
    }

    SmartPtr<Named_interface> source_ni = Root::instance()->interface(
                            gridModels_manager + "/" + grid_name);
    source_grid_ = dynamic_cast<Geostat_grid*>(source_ni.raw_ptr());
    source_rgrid_ = dynamic_cast<RGrid*>(source_ni.raw_ptr());
    Cartesian_grid* source_cartesian = dynamic_cast<Cartesian_grid*>(source_ni.raw_ptr());
    const bool source_is_regular = source_cartesian != NULL;
    if(!source_grid_)
    {
        std::ostringstream error_stream;
        error_stream << grid_name << " is not a valid geostatistical grid object";
        error_mesgs->report("Source_Grid", error_stream.str());
        append_hosc_window_trace("get_Source_Data.fail invalid_geostat_grid");
        return false;
    }
	if(source_cartesian == NULL)
	{
		error_mesgs->report("Source_Grid", "HOSC_WINDOW requires a Cartesian_grid; point sets and irregular grids are not supported");
		append_hosc_window_trace("get_Source_Data.fail non_cartesian_source");
		return false;
	}
    append_hosc_window_tracef(
        "get_Source_Data.grid class=%s is_regular=%d raw_size=%d",
        source_grid_->classname().c_str(),
        source_is_regular ? 1 : 0,
        source_grid_->size());

    string prop_name = parameter_value_compat(parameters, "Source_Property.value");
    append_hosc_window_tracef("get_Source_Data.props requested=%s", prop_name.c_str());
    if(prop_name.empty())
    {
        error_mesgs->report("Source_Property", "No source property specified");
        append_hosc_window_trace("get_Source_Data.fail empty_source_property");
        return false;
    }
    source_val = source_grid_->select_property(prop_name);
    if(source_val == NULL)
    {
        error_mesgs->report("Source_Property", "Selected source property was not found on the source grid");
        append_hosc_window_trace("get_Source_Data.fail missing_source_value_property");
        return false;
    }
    const int input_count = source_grid_->size();
    if(input_count <= 0)
    {
        error_mesgs->report("Source_Grid", "Source grid has no nodes");
        append_hosc_window_trace("get_Source_Data.fail empty_source_grid");
        return false;
    }

	{
		const Geostat_grid::location_type first_location = source_grid_->location(0);
		const GsTLCoordVector source_spacing = source_cartesian->cell_dimensions();
		source_origin_x_ = first_location.x();
		source_origin_y_ = first_location.y();
		source_origin_z_ = first_location.z();
		source_spacing_x_ = source_spacing.x();
		source_spacing_y_ = source_spacing.y();
		source_spacing_z_ = source_spacing.z();
	}
    int min_x = INT_MAX;
    int min_y = INT_MAX;
    int min_z = INT_MAX;
    int max_x = INT_MIN;
    int max_y = INT_MIN;
    int max_z = INT_MIN;
    const int safe_raw_nx = source_cartesian->nx();
    const int safe_raw_ny = source_cartesian->ny();
    const int raw_xy = safe_raw_nx * safe_raw_ny;
    const int safe_raw_nz = source_cartesian->nz();
    if(safe_raw_nx <= 0 || safe_raw_ny <= 0 || safe_raw_nz <= 0 ||
       safe_raw_nx * safe_raw_ny * safe_raw_nz != input_count)
    {
        error_mesgs->report("Source_Grid", "Cartesian grid dimensions do not match its node count");
        append_hosc_window_trace("get_Source_Data.fail inconsistent_cartesian_geometry");
        return false;
    }

    for(int i = 0; i < input_count; i++)
    {
        const int zero_index = i;
        int x = zero_index % safe_raw_nx + 1;
        int y = (zero_index / safe_raw_nx) % safe_raw_ny + 1;
        int z = raw_xy > 0 ? zero_index / raw_xy + 1 : 1;
        if(x < min_x) min_x = x;
        if(y < min_y) min_y = y;
        if(z < min_z) min_z = z;
        if(x > max_x) max_x = x;
        if(y > max_y) max_y = y;
        if(z > max_z) max_z = z;
    }

    if(max_x < min_x || max_y < min_y || max_z < min_z)
    {
        error_mesgs->report("Source_Grid", "Unable to infer source grid dimensions");
        append_hosc_window_trace("get_Source_Data.fail infer_dims_failed");
        return false;
    }

    source_nx_ = max_x - min_x + 1;
    source_ny_ = max_y - min_y + 1;
    source_nz_ = max_z - min_z + 1;
    source_min_x_ = min_x;
    source_min_y_ = min_y;
    source_min_z_ = min_z;
    const int dense_count = source_nx_ * source_ny_ * source_nz_;
    if(source_nx_ <= 0 || source_ny_ <= 0 || source_nz_ <= 0 || dense_count <= 0)
    {
        error_mesgs->report("Source_Grid", "Inferred source grid dimensions are invalid");
        append_hosc_window_tracef("get_Source_Data.fail invalid_dense_dims dims=%d,%d,%d", source_nx_, source_ny_, source_nz_);
        return false;
    }

    append_hosc_window_tracef(
        "get_Source_Data.cartesian dimensions=%d,%d,%d dense_count=%d input_count=%d first_location=%.17g,%.17g,%.17g spacing=%.17g,%.17g,%.17g",
        source_nx_, source_ny_, source_nz_, dense_count, input_count,
        source_origin_x_, source_origin_y_, source_origin_z_,
        source_spacing_x_, source_spacing_y_, source_spacing_z_);

    grid = new (std::nothrow) double*[1 + dense_count];
    if(grid == NULL)
    {
        error_mesgs->report("Source_Grid", "Failed to allocate source grid cache");
        append_hosc_window_trace("get_Source_Data.fail grid_alloc_null");
        return false;
    }
    grid_alloc_rows_ = 0;
    grid[0] = new (std::nothrow) double[2];
    if(grid[0] == NULL)
    {
        error_mesgs->report("Source_Grid", "Failed to allocate source grid header");
        append_hosc_window_trace("get_Source_Data.fail header_alloc_null");
        release_source_grid_cache();
        return false;
    }
    grid_alloc_rows_ = 1;
    grid[0][0] = (double)dense_count;
    grid[0][1] = 4;
    source_present_.assign((size_t)dense_count, 0);

    for(int i = 1; i <= dense_count; i++)
    {
        grid[i] = new (std::nothrow) double[4];
        if(grid[i] == NULL)
        {
            error_mesgs->report("Source_Grid", "Failed to allocate source grid row");
            append_hosc_window_tracef("get_Source_Data.fail row_alloc_null i=%d", i);
            release_source_grid_cache();
            return false;
        }
        grid_alloc_rows_ = i + 1;
        const int zero_index = i - 1;
        const int local_x = zero_index % source_nx_;
        const int local_y = (zero_index / source_nx_) % source_ny_;
        const int local_z = zero_index / (source_nx_ * source_ny_);
        grid[i][0] = (double)(local_x + 1);
        grid[i][1] = (double)(local_y + 1);
        grid[i][2] = (double)(local_z + 1);
        grid[i][3] = 0.0;
    }

    int filled_count = 0;
    for(int i = 0; i < input_count; i++)
    {
		if(!source_val->is_informed(i))
		{
			continue;
		}
        const int zero_index = i;
		const double source_value = (double)source_val->get_value(i);
		if(!wmc::is_finite_number(source_value))
		{
			continue;
		}
        const int dense_index = zero_index + 1;
        grid[dense_index][3] = source_value;
        source_present_[(size_t)(dense_index - 1)] = 1;
        filled_count++;
    }
	if(filled_count == 0)
	{
		error_mesgs->report("Source_Property", "Selected source property has no informed values");
		append_hosc_window_trace("get_Source_Data.fail no_informed_values");
		return false;
	}
	source_grid_name_ = grid_name;
	source_property_name_ = prop_name;
	if((source_nx_ == 1 || source_ny_ == 1) && source_nz_ > 1)
	{
		error_mesgs->report("Source_Grid", "HOSC_WINDOW currently supports XY single-layer or full 3D Cartesian grids; XZ/YZ degenerate grids are not supported");
		append_hosc_window_tracef("get_Source_Data.fail unsupported_degenerate_axes dims=%d,%d,%d", source_nx_, source_ny_, source_nz_);
		return false;
	}

    append_hosc_window_tracef(
        "get_Source_Data.end dense_rows=%d filled=%d missing_zero=%d first=(%.17g,%.17g,%.17g,%.17g) last=(%.17g,%.17g,%.17g,%.17g)",
        dense_count,
        filled_count,
        dense_count - filled_count,
        grid[1][0], grid[1][1], grid[1][2], grid[1][3],
        grid[dense_count][0], grid[dense_count][1], grid[dense_count][2], grid[dense_count][3]);
    return true;
}

bool HOSC_WINDOW::get_Target_Data(const Parameters_handler* parameters,
							         Error_messages_handler* error_mesgs )
{

    string grid_name = parameter_value_compat(parameters, "Target_Grid.value" );

    if ( grid_name.empty() )
    {
        error_mesgs->report( "Target_Grid", "No rotation grid selected" );
        return false;
    }


    target_grid_ = dynamic_cast<RGrid*>( Root::instance()->interface(
                            gridModels_manager + "/" + grid_name).raw_ptr() );

    if( !target_grid_ )
    {
        std::ostringstream error_stream;
        error_stream <<  grid_name <<  " is not a valid (Regular) object grid";
        error_mesgs->report( "Target_Grid", error_stream.str() );
        return false;
    }











	target_prop_name_ = parameter_value_compat(parameters, "props.value" );
	if(target_prop_name_.empty())
	{
		error_mesgs->report("props","No property specified");
		return false;
	}

#if 0
	target_dim_num = parameter_value_compat(parameters, "dims.value");
	if(target_dim_num.empty())
	{
		error_mesgs->report("dims", "No dimmension specified");
		return false;
	}
	const char* tempDimStr = target_dim_num.c_str();
	int dimnum = atoi(tempDimStr);
	if(dimnum > Key->NDIR || dimnum < 2 || dimnum > 3)
	{
		error_mesgs->report("dims", "Wrong dimension");
		return false;
	}
#endif
	target_property_ = geostat_utils::add_property_to_grid(target_grid_, target_prop_name_);

    return true;
}

bool HOSC_WINDOW::get_Keys_Data(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	Key = (Keys*)calloc(1, sizeof(Keys));








	Key->REG_IRREG = 1;




	Key->NX = source_nx_;




	Key->NY = source_ny_;




	Key->NZ = source_nz_;

	if(Key->NZ == 0)
		Key->DIM = 2;
	else
		Key->DIM = 3;

#if 0
	string DIMX_TEXT = parameter_value_compat(parameters, "DIMX.value");
	stringstream ss6(DIMX_TEXT);
	ss6 >> Key->DIMX;

	string DIMY_TEXT = parameter_value_compat(parameters, "DIMY.value");
	stringstream ss7(DIMY_TEXT);
	ss7 >> Key->DIMY;

	string DIMZ_TEXT = parameter_value_compat(parameters, "DIMZ.value");
	stringstream ss8(DIMZ_TEXT);
	ss8 >> Key->DIMZ;
#endif
	if(Key->NX > 0)
		Key->DIMX = Key->NX-1;
	else
		Key->DIMX = 0;
	if(Key->NY > 0)
		Key->DIMY = Key->NY-1;
	else
		Key->DIMY = 0;
	if(Key->NZ > 0)
		Key->DIMZ = Key->NZ-1;
	else
		Key->DIMZ = 0;
	MaxX = Key->NX;
	MaxY = Key->NY;
	MaxZ = Key->NZ;


	string ORDER_TEXT = parameter_value_compat(parameters, "ORDER.value");
	stringstream ss9(ORDER_TEXT);
	ss9 >> Key->ORDER;




	Key->NDIR = Key->ORDER-1;








	vector<string>texts;
	texts.push_back(parameter_value_compat(parameters, "NOD.value"));
	texts.push_back(parameter_value_compat(parameters, "NLAG.value"));
	texts.push_back(parameter_value_compat(parameters, "LAG.value"));
	texts.push_back(parameter_value_compat(parameters, "LTOL.value"));
	texts.push_back(parameter_value_compat(parameters, "AZM.value"));
	texts.push_back(parameter_value_compat(parameters, "ATOL.value"));
	texts.push_back(parameter_value_compat(parameters, "BANDWH.value"));
	texts.push_back(parameter_value_compat(parameters, "DIP.value"));
	texts.push_back(parameter_value_compat(parameters, "DTOL.value"));
	texts.push_back(parameter_value_compat(parameters, "BANDWD.value"));


	vector<string> tempResult;
	for(int i = 0; i < texts.size(); i++)
	{
		stringstream ss(texts[i]);
		string strItem;
		int tempNum = 0;
		while(getline(ss, strItem, ','))
		{
			tempResult.push_back(strItem);
			tempNum++;
		}
		if(tempNum != Key->NDIR)
		{

			error_mesgs->report("NLAG", "Number of NLAG should equal to number of direction");
			return false;
		}
	}

	vector<int*>tempInt;
	Key->NOD = (int*)calloc(Key->NDIR, sizeof(int));
	tempInt.push_back(Key->NOD);
	Key->NLAG = (int*)calloc(Key->NDIR, sizeof(int));
	tempInt.push_back(Key->NLAG);
	Key->LAG = (int*)calloc(Key->NDIR, sizeof(int));
	tempInt.push_back(Key->LAG);
	Key->LTOL = (int*)calloc(Key->NDIR, sizeof(int));
	tempInt.push_back(Key->LTOL);
	Key->AZM = (int*)calloc(Key->NDIR, sizeof(int));
	tempInt.push_back(Key->AZM);
	Key->ATOL = (int*)calloc(Key->NDIR, sizeof(int));
	tempInt.push_back(Key->ATOL);
	Key->BANDWH = (int*)calloc(Key->NDIR, sizeof(int));
	tempInt.push_back(Key->BANDWH);
	Key->DIP = (int*)calloc(Key->NDIR, sizeof(int));
	tempInt.push_back(Key->DIP);
	Key->DTOL = (int*)calloc(Key->NDIR, sizeof(int));
	tempInt.push_back(Key->DTOL);
	Key->BANDWD = (int*)calloc(Key->NDIR, sizeof(int));
	tempInt.push_back(Key->BANDWD);

	for(int i = 0; i < tempResult.size(); i++)
	{
		for(int j = 0; j < Key->NDIR; j++)
		{
			tempInt[i/Key->NDIR][j] = atoi(tempResult[i].c_str());
			i++;
		}
		i--;
	}
	vector<string>().swap(texts);
	vector<string>().swap(tempResult);







	return true;
}


int HOSC_WINDOW::set_Target_Data(double** targetGrid, int gridSize, int valuePos)
{

	for(int i = 0; i < gridSize; i++)
	{
		target_property_->set_value(targetGrid[i][valuePos], i);

	}


	return 0;
}


bool HOSC_WINDOW::get_Calculate_Target(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	string Calculate_Target = parameter_value_compat(parameters, "Target_Button.value" );
	append_hosc_window_tracef("get_Calculate_Target raw=%s", Calculate_Target.c_str());

	if(Calculate_Target == "Moment")
		G_CalculateTarget = 3;
	else if(Calculate_Target == "Cumulant")
		G_CalculateTarget = 4;
	else
	{
		error_mesgs->report("Target_Button", "Expected Moment or Cumulant");
		append_hosc_window_trace("get_Calculate_Target.fail unknown_target");
		return false;
	}

	return true;
}

bool HOSC_WINDOW::get_Calculate_Algorithm(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	append_hosc_window_trace("get_Calculate_Algorithm.before_value");
	string Calculate_Algorithm = parameter_value_compat(parameters, "Algorithm_Button.value");
	append_hosc_window_trace("get_Calculate_Algorithm.enter");
	append_hosc_window_tracef("get_Calculate_Algorithm.raw_length=%d", (int)Calculate_Algorithm.size());
	if(Calculate_Algorithm == "Bell" || Calculate_Algorithm == "bell" || Calculate_Algorithm.find("Bell") != string::npos || Calculate_Algorithm.find("bell") != string::npos)
		al = Bell;
	else if(Calculate_Algorithm == "AVL" || Calculate_Algorithm == "avl" || Calculate_Algorithm.find("AVL") != string::npos || Calculate_Algorithm.find("avl") != string::npos)
		al = AVL;
	else
	{
		error_mesgs->report("Algorithm_Button", "Expected Bell or AVL");
		append_hosc_window_trace("get_Calculate_Algorithm.fail unknown_algorithm");
		return false;
	}
	append_hosc_window_tracef("get_Calculate_Algorithm.parsed algorithm=%d", (int)al);
	append_hosc_window_trace("get_Calculate_Algorithm.before_return");
	return true;
}

bool HOSC_WINDOW::get_Weighted_Moment_Options(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	append_hosc_window_trace("get_Weighted_Moment_Options.begin");
	weightedOptions = wmc::WeightedMomentOptions();
	weightedOptions.normalize_weights = true;
	weightedOptions.zero_first_order_moment = false;
	useWeightedMoment = false;
	coordinateWeightUsesWindowCenter = true;

	string enabled = parameter_value_compat(parameters, "enable_weighted_moment.value");
	if(enabled.empty())
	{
		enabled = parameter_value_compat(parameters, "Use_Weighted_Moment.value");
	}
	if(!enabled.empty() && !parse_strict_bool_text(enabled, &useWeightedMoment))
	{
		error_mesgs->report("enable_weighted_moment", "Expected a boolean value");
		return false;
	}

	const string mode = parameter_value_compat(parameters, "weight_mode.value");
	if(!mode.empty() && !parse_weight_mode_text(mode, &weightedOptions.mode))
	{
		error_mesgs->report("weight_mode", "Supported modes are uniform, coordinate_gaussian, and soft_high_porosity");
		return false;
	}
	if(!useWeightedMoment && weightedOptions.mode != wmc::WeightMode_Uniform)
	{
		error_mesgs->report("enable_weighted_moment", "A non-uniform weight mode requires weighted calculation to be enabled");
		return false;
	}

	const string bandwidth = parameter_value_compat(parameters, "weight_bandwidth.value");
	if(!bandwidth.empty() && !parse_strict_double_text(bandwidth, &weightedOptions.bandwidth))
	{
		error_mesgs->report("weight_bandwidth", "Expected a finite number");
		return false;
	}
	if(weightedOptions.mode == wmc::WeightMode_CoordinateGaussian &&
	   weightedOptions.bandwidth <= 0.0)
	{
		error_mesgs->report("weight_bandwidth", "Coordinate Gaussian sigma must be greater than zero");
		return false;
	}

	const string min_weight = parameter_value_compat(parameters, "min_weight.value");
	if(!min_weight.empty() && !parse_strict_double_text(min_weight, &weightedOptions.min_weight))
	{
		error_mesgs->report("min_weight", "Expected a finite number");
		return false;
	}
	if(weightedOptions.min_weight < 0.0)
	{
		error_mesgs->report("min_weight", "Minimum weight must not be negative");
		return false;
	}

	const string normalize = parameter_value_compat(parameters, "normalize_weights.value");
	bool requested_normalize = true;
	if(!normalize.empty() && !parse_strict_bool_text(normalize, &requested_normalize))
	{
		error_mesgs->report("normalize_weights", "Expected a boolean value");
		return false;
	}
	if(!requested_normalize)
	{
		error_mesgs->report("normalize_weights", "Window weights are always normalized. Remove normalize_weights from the parameter file or set it to true.");
		return false;
	}

	const string zero_first = parameter_value_compat(parameters, "zero_first_order_moment.value");
	bool requested_zero_first = false;
	if(!zero_first.empty() && !parse_strict_bool_text(zero_first, &requested_zero_first))
	{
		error_mesgs->report("zero_first_order_moment", "Expected a boolean value");
		return false;
	}
	if(requested_zero_first)
	{
		error_mesgs->report("zero_first_order_moment", "Legacy zero-first-order centering is disabled. Remove zero_first_order_moment from the parameter file or set it to false; use center_cumulant_input_for_order_ge_2 for cumulants.");
		return false;
	}

	const string center_mode = parameter_value_compat(parameters, "weight_center_mode.value");
	if(!center_mode.empty())
	{
		const string normalized_center_mode = normalize_algorithm_text(center_mode);
		if(!normalized_center_mode.empty() && normalized_center_mode != "window_center" &&
		   normalized_center_mode != "window center")
		{
			error_mesgs->report("weight_center_mode", "Coordinate Gaussian weights always use the window center. Remove weight_center_mode from the parameter file or set it to window_center.");
			append_hosc_window_trace("get_Weighted_Moment_Options.fail parameter_center_rejected");
			return false;
		}
	}

	if(weightedOptions.mode == wmc::WeightMode_SoftHighPorosity)
	{
		const string threshold = parameter_value_compat(parameters, "high_porosity_threshold.value");
		const string transition = parameter_value_compat(parameters, "soft_transition_width.value");
		const string strength = parameter_value_compat(parameters, "soft_weight_strength.value");
		if(!threshold.empty() && !parse_strict_double_text(threshold, &weightedOptions.high_porosity_threshold))
		{
			error_mesgs->report("high_porosity_threshold", "Expected a finite number");
			return false;
		}
		if(!transition.empty() && !parse_strict_double_text(transition, &weightedOptions.soft_transition_width))
		{
			error_mesgs->report("soft_transition_width", "Expected a finite number");
			return false;
		}
		if(!strength.empty() && !parse_strict_double_text(strength, &weightedOptions.soft_weight_strength))
		{
			error_mesgs->report("soft_weight_strength", "Expected a finite number");
			return false;
		}
		if(weightedOptions.soft_transition_width <= 0.0 ||
		   weightedOptions.soft_weight_strength < 0.0 ||
		   weightedOptions.soft_weight_strength > 1.0)
		{
			error_mesgs->report("soft_weight_strength", "Soft weighting requires transition width > 0 and strength in [0, 1]");
			return false;
		}
	}
	append_hosc_window_tracef(
		"get_Weighted_Moment_Options.end enabled=%d mode=%d window_center=%d bandwidth=%.17g threshold=%.17g transition=%.17g strength=%.17g min=%.17g normalize=%d zero_first=%d",
		useWeightedMoment ? 1 : 0,
		(int)weightedOptions.mode,
		coordinateWeightUsesWindowCenter ? 1 : 0,
		weightedOptions.bandwidth,
		weightedOptions.high_porosity_threshold,
		weightedOptions.soft_transition_width,
		weightedOptions.soft_weight_strength,
		weightedOptions.min_weight,
		weightedOptions.normalize_weights ? 1 : 0,
		weightedOptions.zero_first_order_moment ? 1 : 0);
	return true;
}

#if 0
bool HOSC::get_Dim_Data(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	string Dim = parameter_value_compat(parameters, "dims.value");
	const char* cDim = Dim.c_str();
	if(Dim.empty())
	{
		error_mesgs->report("dims", "input dimension of result");
		G_TargetDimension = 0;
		return false;
	}
	int dimension = atoi(cDim);
	if(dimension > 3 || dimension < 1)
	{
		error_mesgs->report("dims", "dimension out of range (1~3)");
		G_TargetDimension = 0;
		return false;
	}
	G_TargetDimension = dimension;
	return true;
}
#endif





bool HOSC_WINDOW::get_Windows_Moment(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{

	append_hosc_window_trace("get_Windows_Moment.begin");

	string wl = parameter_value_compat(parameters, "windowl.value");
	string ww = parameter_value_compat(parameters, "windoww.value");
	string wh = parameter_value_compat(parameters, "windowh.value");

	string order = parameter_value_compat(parameters, "moment_order.value");

	string sx = parameter_value_compat(parameters, "stepx.value");
	string sy = parameter_value_compat(parameters, "stepy.value");
	string sz = parameter_value_compat(parameters, "stepz.value");
	string invert_output = parameter_value_compat(parameters, "invert_window_output.value");
	string enable_clip = parameter_value_compat(parameters, "enable_percentile_clipping.value");
	string clip_min = parameter_value_compat(parameters, "percentile_vmin.value");
	string clip_max = parameter_value_compat(parameters, "percentile_vmax.value");
	string diagnostics = parameter_value_compat(parameters, "output_diagnostics.value");
	string center_cumulant = parameter_value_compat(parameters, "center_cumulant_input_for_order_ge_2.value");
	string standardized = parameter_value_compat(parameters, "output_standardized_cumulant.value");
	string variance_epsilon = parameter_value_compat(parameters, "variance_epsilon.value");
	string write_audit_sidecar = parameter_value_compat(parameters, "write_audit_sidecar.value");
	string audit_output_directory = parameter_value_compat(parameters, "audit_output_directory.value");
	invertWindowOutput = false;
	enablePercentileClipping = false;
	outputDiagnostics = true;
	centerCumulantInputForOrderGE2 = true;
	outputStandardizedCumulant = false;
	writeAuditSidecar = false;
	percentileClipMin = 2.0;
	percentileClipMax = 98.0;
	varianceEpsilon = 1.0e-12;
	auditOutputDirectory.clear();
	auditRunDirectory.clear();
	bool requested_invert = false;
	if(!invert_output.empty() && !parse_strict_bool_text(invert_output, &requested_invert))
	{
		error_mesgs->report("invert_window_output", "Expected a boolean value");
		return false;
	}
	if(requested_invert)
	{
		error_mesgs->report("invert_window_output", "The primary statistic is never inverted. Remove invert_window_output from the parameter file or set it to false.");
		return false;
	}
	if(!enable_clip.empty() && !parse_strict_bool_text(enable_clip, &enablePercentileClipping))
	{
		error_mesgs->report("enable_percentile_clipping", "Expected a boolean value");
		return false;
	}
	if(!diagnostics.empty() && !parse_strict_bool_text(diagnostics, &outputDiagnostics))
	{
		error_mesgs->report("output_diagnostics", "Expected a boolean value");
		return false;
	}
	if(!center_cumulant.empty() && !parse_strict_bool_text(center_cumulant, &centerCumulantInputForOrderGE2))
	{
		error_mesgs->report("center_cumulant_input_for_order_ge_2", "Expected a boolean value");
		return false;
	}
	if(!standardized.empty() && !parse_strict_bool_text(standardized, &outputStandardizedCumulant))
	{
		error_mesgs->report("output_standardized_cumulant", "Expected a boolean value");
		return false;
	}
	if(!write_audit_sidecar.empty() && !parse_strict_bool_text(write_audit_sidecar, &writeAuditSidecar))
	{
		error_mesgs->report("write_audit_sidecar", "Expected a boolean value");
		return false;
	}
	auditOutputDirectory = audit_output_directory;
	auditRunDirectory.clear();
	if(writeAuditSidecar && auditOutputDirectory.empty())
	{
		error_mesgs->report("audit_output_directory", "An audit output directory is required when write_audit_sidecar is true");
		return false;
	}
	if(!variance_epsilon.empty() && !parse_strict_double_text(variance_epsilon, &varianceEpsilon))
	{
		error_mesgs->report("variance_epsilon", "Expected a finite number");
		return false;
	}
	if(varianceEpsilon < 0.0)
	{
		error_mesgs->report("variance_epsilon", "Variance epsilon must not be negative");
		return false;
	}
	if(!clip_min.empty() && !parse_strict_double_text(clip_min, &percentileClipMin))
	{
		error_mesgs->report("percentile_vmin", "Expected a finite number");
		return false;
	}
	if(!clip_max.empty() && !parse_strict_double_text(clip_max, &percentileClipMax))
	{
		error_mesgs->report("percentile_vmax", "Expected a finite number");
		return false;
	}
	append_hosc_window_tracef(
		"get_Windows_Moment.raw window=%s,%s,%s order=%s step=%s,%s,%s invert=%s clip=%s clip_range=%s,%s",
		wl.c_str(), ww.c_str(), wh.c_str(), order.c_str(), sx.c_str(), sy.c_str(), sz.c_str(), invert_output.c_str(), enable_clip.c_str(), clip_min.c_str(), clip_max.c_str());




	auto_prop_name_ = parameter_value_compat(parameters, "props.value" );
	if(auto_prop_name_.empty())
	{
		error_mesgs->report("props","No property specified");
		append_hosc_window_trace("get_Windows_Moment.fail empty_target_property");
		return false;
	}

	int iwl = 0;
	int iww = 0;
	int iwh = 0;
	int io = 0;
	if(!parse_strict_int_text(wl, &iwl) || !parse_strict_int_text(ww, &iww) ||
	   !parse_strict_int_text(wh, &iwh) || !parse_strict_int_text(order, &io) ||
	   !parse_strict_int_text(sx, &stepX) || !parse_strict_int_text(sy, &stepY) ||
	   !parse_strict_int_text(sz, &stepZ))
	{
		error_mesgs->report("moment_order", "Window sizes, order, and move steps must be complete integers");
		append_hosc_window_trace("get_Windows_Moment.fail non_integer_parameter");
		return false;
	}

	if(iwl <= 0 || iww <= 0)
	{
		error_mesgs->report("windowl","Window X and Y must be positive");
		append_hosc_window_tracef("get_Windows_Moment.fail invalid_window_xy wx=%d wy=%d", iwl, iww);
		return false;
	}
	if(stepX <= 0 || stepY <= 0)
	{
		error_mesgs->report("stepx","Move Step X and Y must be positive");
		append_hosc_window_tracef("get_Windows_Moment.fail invalid_step_xy sx=%d sy=%d", stepX, stepY);
		return false;
	}
	if(source_nx_ > 1 && source_ny_ > 1 && source_nz_ > 1 && stepZ <= 0)
	{
		error_mesgs->report("stepz", "Move Step Z must be positive for 3D grids");
		append_hosc_window_tracef("get_Windows_Moment.fail invalid_step_z sz=%d", stepZ);
		return false;
	}
	if(source_nx_ > 1 && source_ny_ > 1 && source_nz_ > 1 && iwh <= 0)
	{
		error_mesgs->report("windowh", "Window Z must be positive for 3D grids");
		append_hosc_window_tracef("get_Windows_Moment.fail invalid_window_z wz=%d", iwh);
		return false;
	}


	if(source_nz_ == 1)
	{
		if((iwh != 1 && iwh != 0))
		{
			error_mesgs->report("windowz","should be 1 if souce grid is 2D");
			append_hosc_window_tracef("get_Windows_Moment.fail invalid_2d_window_z wz=%d", iwh);
			return false;
		}
		if(stepZ != 0)
		{
			error_mesgs->report("stepz","should be 0 if souce grid is 2D");
			append_hosc_window_tracef("get_Windows_Moment.fail invalid_2d_step_z sz=%d", stepZ);
			return false;
		}
		iwh = 1;
	}


	momentWindowL = iwl;
	momentWindowW = iww;
	momentWindowH = iwh;
	momentOrder = io;
	if(percentileClipMin > percentileClipMax)
	{
		const double temp_percentile = percentileClipMin;
		percentileClipMin = percentileClipMax;
		percentileClipMax = temp_percentile;
	}
	if(momentOrder < 1 || momentOrder > 8)
	{
		error_mesgs->report("moment_order","Supported orders are integers from 1 through 8");
		append_hosc_window_tracef("get_Windows_Moment.fail invalid_order order=%d", momentOrder);
		return false;
	}
	if(outputStandardizedCumulant && G_CalculateTarget != 4)
	{
		error_mesgs->report("output_standardized_cumulant", "Standardized cumulant output requires Target_Button=Cumulant");
		return false;
	}

	int iinx = 0;
	int iiny = 0;
	int iinz = 0;

	if(momentWindowL > source_nx_ || momentWindowW > source_ny_ ||
	   (source_nz_ > 1 && momentWindowH > source_nz_))
	{
		error_mesgs->report("windowl", "A complete moving window must fit inside the source grid");
		append_hosc_window_trace("get_Windows_Moment.fail full_window_does_not_fit");
		return false;
	}

	iinx = 1 + (source_nx_ - momentWindowL) / stepX;
	iiny = 1 + (source_ny_ - momentWindowW) / stepY;
	if(source_nx_ == 1 || source_ny_ == 1 || source_nz_ == 1)
		iinz = 1;
	else
		iinz = 1 + (source_nz_ - momentWindowH) / stepZ;

	if(iinx <= 0 || iiny <= 0)
	{
		error_mesgs->report("windowl","No less than 0");
		append_hosc_window_tracef("get_Windows_Moment.fail invalid_output_xy out=%d,%d", iinx, iiny);
		return false;
	}
	if(source_nx_ > 1 && source_ny_ > 1 && source_nz_ > 1 && iinz <= 0)
	{
		error_mesgs->report("windowH", "No less than 0");
	}

	double output_dz = (source_nx_ == 1 || source_ny_ == 1 || source_nz_ == 1) ? 1.0 : (double)stepZ;
	const double center_x_origin = source_origin_x_ + 0.5 * (double)(momentWindowL - 1) * source_spacing_x_;
	const double center_y_origin = source_origin_y_ + 0.5 * (double)(momentWindowW - 1) * source_spacing_y_;
	const double center_z_origin = source_origin_z_ + 0.5 * (double)(momentWindowH - 1) * source_spacing_z_;
	GridShape tempShape = {
		iinx, iiny, iinz,
		center_x_origin, center_y_origin, center_z_origin,
		(double)stepX * source_spacing_x_,
		(double)stepY * source_spacing_y_,
		output_dz * source_spacing_z_};


	myWindowShape = tempShape;
	if(writeAuditSidecar)
	{
		if(!create_window_audit_directory(
			auditOutputDirectory,
			processId_,
			runSequence_,
			&auditRunDirectory))
		{
			error_mesgs->report("audit_output_directory", "Unable to create a unique audit run directory");
			append_hosc_window_instance_trace(this, runSequence_, "audit.start_directory_failed");
			return false;
		}
		if(!write_window_audit_status(
			auditRunDirectory,
			runId_,
			processId_,
			runSequence_,
			"started",
			"parameters_validated"))
		{
			error_mesgs->report("audit_output_directory", "Unable to write the audit run status");
			append_hosc_window_instance_trace(this, runSequence_, "audit.start_status_failed");
			return false;
		}
		append_hosc_window_tracef(
			"run_sequence=%lu instance=%p stage=audit.started directory=%s",
			runSequence_,
			this,
			auditRunDirectory.c_str());
	}
	append_hosc_window_tracef(
		"get_Windows_Moment.end target_prop=%s parsed_window=%d,%d,%d parsed_step=%d,%d,%d output=%d,%d,%d invert=%d clip=%d diagnostics=%d center_cumulant=%d standardized=%d variance_epsilon=%.17g clip_percentiles=%.17g,%.17g",
		auto_prop_name_.c_str(), momentWindowL, momentWindowW, momentWindowH, stepX, stepY, stepZ, myWindowShape.nx, myWindowShape.ny, myWindowShape.nz, invertWindowOutput ? 1 : 0, enablePercentileClipping ? 1 : 0, outputDiagnostics ? 1 : 0, centerCumulantInputForOrderGE2 ? 1 : 0, outputStandardizedCumulant ? 1 : 0, varianceEpsilon, percentileClipMin, percentileClipMax);
	append_hosc_window_instance_trace(this, runSequence_, "window_params.body.end.locals_destroy_next");
	return true;
}


int HOSC_WINDOW::calculateMomentWindow()
{
	if(momentWindowL == 0 && momentWindowW == 0 && momentWindowH == 0)
	{
		append_hosc_window_trace("calculateMomentWindow: all window dimensions are zero");
		return -1;
	}

	const int source_nx = source_nx_;
	const int source_ny = source_ny_;
	const int source_nz = source_nz_;
	const bool is_2d = is_effective_2d_dims(source_nx, source_ny, source_nz);
	const int out_nx = myWindowShape.nx;
	const int out_ny = myWindowShape.ny;
	const int out_nz = myWindowShape.nz;
	const int gridSize = out_nx * out_ny * out_nz;
	const int output_effective_dim =
		((out_nx > 1) ? 1 : 0) +
		((out_ny > 1) ? 1 : 0) +
		((out_nz > 1) ? 1 : 0);
	G_TargetDimension = is_2d ? 2 : 3;
	if(gridSize <= 0)
	{
		append_hosc_window_trace("calculateMomentWindow: invalid output grid size");
		return -1;
	}

	append_hosc_window_tracef(
		"calculateMomentWindow.begin build=movewindow_full_window_cartesian_20260726 target=%d algorithm=%d weighted=%d weight_mode=%d invert=%d clip=%d clip_percentiles=%.17g,%.17g order=%d window=%d,%d,%d step=%d,%d,%d source=%d,%d,%d output=%d,%d,%d output_dim=%d",
		G_CalculateTarget,
		(int)al,
		useWeightedMoment ? 1 : 0,
		(int)weightedOptions.mode,
		invertWindowOutput ? 1 : 0,
		enablePercentileClipping ? 1 : 0,
		percentileClipMin,
		percentileClipMax,
		momentOrder,
		momentWindowL,
		momentWindowW,
		momentWindowH,
		stepX,
		stepY,
		stepZ,
		source_nx,
		source_ny,
		source_nz,
		out_nx,
		out_ny,
		out_nz,
		output_effective_dim);

	int block = 0;
	int kept_count = 0;
	std::vector<double> output_values;
	std::vector<double> output_support;
	std::vector<double> output_ess;
	std::vector<double> output_numeric_valid;
	std::vector<double> output_display_valid;
	std::vector<double> output_kappa2;
	std::vector<double> output_lambda;
	std::vector<double> output_lambda_valid;
	output_values.assign((size_t)gridSize, 0.0);
	output_support.assign((size_t)gridSize, 0.0);
	output_ess.assign((size_t)gridSize, 0.0);
	output_numeric_valid.assign((size_t)gridSize, 0.0);
	output_display_valid.assign((size_t)gridSize, 0.0);
	output_kappa2.assign((size_t)gridSize, 0.0);
	output_lambda.assign((size_t)gridSize, 0.0);
	output_lambda_valid.assign((size_t)gridSize, 0.0);
	output_keep_mask_.assign((size_t)gridSize, 0);
	for(int anchor_z = 0; anchor_z < out_nz; anchor_z++)
	{
		const int start_z = is_2d ? 0 : anchor_z * stepZ;
		for(int anchor_y = 0; anchor_y < out_ny; anchor_y++)
		{
			const int start_y = anchor_y * stepY;
			for(int anchor_x = 0; anchor_x < out_nx; anchor_x++)
			{
				const int start_x = anchor_x * stepX;
				const int length = momentWindowL;
				const int width = momentWindowW;
				const int height = is_2d ? 1 : momentWindowH;
				wmc::WeightedMomentOptions window_weight_options = weightedOptions;
				if(useWeightedMoment &&
				   weightedOptions.mode == wmc::WeightMode_CoordinateGaussian &&
				   coordinateWeightUsesWindowCenter)
				{
					// Coordinate Gaussian is defined in relative window coordinates.
					// Each window therefore has its own geometric center and half-width.
					if(!wmc::configure_relative_window_coordinate_gaussian(
						&window_weight_options,
						source_origin_x_, source_origin_y_, source_origin_z_,
						source_spacing_x_, source_spacing_y_, source_spacing_z_,
						start_x, start_y, start_z, length, width, height))
					{
						append_hosc_window_tracef("calculateMomentWindow.invalid_coordinate_geometry block=%d", block);
						return -1;
					}
				}
				if(length <= 0 || width <= 0 || height <= 0 || block >= gridSize)
				{
					append_hosc_window_tracef(
						"calculateMomentWindow.invalid_window block=%d start=%d,%d,%d size=%d,%d,%d",
						block, start_x, start_y, start_z, length, width, height);
					return -1;
				}

				std::vector<double> samples;
				std::vector<double> weights;
				samples.reserve((size_t)length * (size_t)width * (size_t)height);
				if(useWeightedMoment)
				{
					weights.reserve((size_t)length * (size_t)width * (size_t)height);
				}
				bool complete_data_window = true;
				for(int z = 0; z < height; z++)
				{
					for(int y = 0; y < width; y++)
					{
						for(int x = 0; x < length; x++)
						{
							const int gx = start_x + x;
							const int gy = start_y + y;
							const int gz = start_z + z;
							const int source_index = gz * source_nx * source_ny + gy * source_nx + gx + 1;
							if(source_index < 1 || source_index > (int)grid[0][0])
							{
								append_hosc_window_tracef(
									"calculateMomentWindow.source_index_oob block=%d index=%d source_total=%.0f",
									block, source_index, grid[0][0]);
								return -1;
							}
							if(source_index - 1 >= (int)source_present_.size() ||
							   source_present_[(size_t)(source_index - 1)] == 0)
							{
								complete_data_window = false;
								continue;
							}
							const double sample_value = grid[source_index][3];
							samples.push_back(sample_value);
							if(useWeightedMoment)
							{
								wmc::WeightedSample weighted_sample;
								weighted_sample.values = &sample_value;
								weighted_sample.value_count = 1;
								weighted_sample.x = source_origin_x_ + (double)gx * source_spacing_x_;
								weighted_sample.y = source_origin_y_ + (double)gy * source_spacing_y_;
								weighted_sample.z = source_origin_z_ + (double)gz * source_spacing_z_;
								weights.push_back(wmc::compute_weight(weighted_sample, window_weight_options));
							}
						}
					}
				}
				if(!complete_data_window || samples.empty())
				{
					output_support[(size_t)block] = (double)samples.size();
					output_values[(size_t)block] = 0.0;
					block++;
					continue;
				}
				output_support[(size_t)block] = (double)samples.size();
				if(useWeightedMoment)
				{
					const wmc::WeightDiagnostics weight_diagnostics =
						wmc::compute_weight_diagnostics(&weights[0], (int)weights.size());
					output_ess[(size_t)block] = weight_diagnostics.effective_sample_size;
					if(weight_diagnostics.positive_weight_count <= 0 ||
					   !wmc::is_finite_number(weight_diagnostics.weight_sum) ||
					   !wmc::is_finite_number(weight_diagnostics.effective_sample_size))
					{
						append_hosc_window_tracef("calculateMomentWindow.invalid_weights block=%d", block);
						output_values[(size_t)block] = 0.0;
						block++;
						continue;
					}
				}
				else
				{
					output_ess[(size_t)block] = (double)samples.size();
				}

				const bool cumulant_centering =
					G_CalculateTarget == 4 && momentOrder >= 2 && centerCumulantInputForOrderGE2;
				std::vector<double> centered_samples;
				const std::vector<double>* unweighted_cumulant_samples = &samples;
				if(!useWeightedMoment && cumulant_centering)
				{
					centered_samples = centered_unweighted_samples(samples);
					unweighted_cumulant_samples = &centered_samples;
				}
				double value = 0.0;
				double raw_cumulant = 0.0;
				double kappa2 = 0.0;
				bool lambda_valid = false;
				if(useWeightedMoment)
				{
					try
					{
						if(G_CalculateTarget == 3)
						{
							value = wmc::weighted_univariate_moment(
								&samples[0], &weights[0], (int)samples.size(), momentOrder,
								window_weight_options.normalize_weights, window_weight_options.zero_first_order_moment);
						}
						else if(G_CalculateTarget == 4)
						{
							raw_cumulant = wmc::weighted_univariate_cumulant(
								&samples[0], &weights[0], (int)samples.size(), momentOrder,
								window_weight_options.normalize_weights, cumulant_centering);
							value = raw_cumulant;
						}
						else
						{
							append_hosc_window_tracef("calculateMomentWindow.unsupported_weighted_target target=%d", G_CalculateTarget);
							return -1;
						}
					}
					catch(const std::exception& e)
					{
						append_hosc_window_tracef("calculateMomentWindow.weighted_exception block=%d message=%s", block, e.what());
						return -1;
					}
				}
				else if(G_CalculateTarget == 3)
				{
					value = calculate_window_moment_from_samples(samples, momentOrder);
				}
				else if(G_CalculateTarget == 4)
				{
					if(!calculate_window_cumulant_from_samples(*unweighted_cumulant_samples, momentOrder, al, &raw_cumulant))
					{
						append_hosc_window_tracef("calculateMomentWindow.unweighted_cumulant_failed block=%d sample_count=%d", block, (int)samples.size());
						return -1;
					}
					value = raw_cumulant;
				}
				else
				{
					append_hosc_window_tracef("calculateMomentWindow.unsupported_target target=%d", G_CalculateTarget);
					return -1;
				}

				if(outputStandardizedCumulant && G_CalculateTarget == 4)
				{
					try
					{
						if(useWeightedMoment)
						{
							kappa2 = wmc::weighted_univariate_cumulant(
								&samples[0], &weights[0], (int)samples.size(), 2,
								window_weight_options.normalize_weights, cumulant_centering);
						}
						else if(!calculate_window_cumulant_from_samples(*unweighted_cumulant_samples, 2, al, &kappa2))
						{
							append_hosc_window_tracef("calculateMomentWindow.kappa2_failed block=%d", block);
							return -1;
						}
					}
					catch(const std::exception& e)
					{
						append_hosc_window_tracef("calculateMomentWindow.kappa2_exception block=%d message=%s", block, e.what());
						return -1;
					}
					output_kappa2[(size_t)block] = kappa2;
					if(wmc::is_finite_number(raw_cumulant) && wmc::is_finite_number(kappa2) && kappa2 > varianceEpsilon)
					{
						const double lambda = raw_cumulant / pow(kappa2, 0.5 * (double)momentOrder);
						if(wmc::is_finite_number(lambda))
						{
							output_lambda[(size_t)block] = lambda;
							output_lambda_valid[(size_t)block] = 1.0;
							lambda_valid = true;
						}
					}
				}

				if(!wmc::is_finite_number(value))
				{
					append_hosc_window_tracef("calculateMomentWindow.nonfinite_result block=%d", block);
					output_values[(size_t)block] = 0.0;
					block++;
					continue;
				}

				output_keep_mask_[(size_t)block] = 1;
				output_numeric_valid[(size_t)block] = 1.0;
				if(is_float_representable(value))
				{
					output_display_valid[(size_t)block] = 1.0;
				}
				else
				{
					append_hosc_window_tracef("calculateMomentWindow.float_overflow block=%d value=%.17g", block, value);
				}
				kept_count++;

				if(block == 0)
				{
					append_hosc_window_tracef(
						"calculateMomentWindow.first_block start=%d,%d,%d size=%d,%d,%d samples=%d weights=%d ess=%.17g value=%.17g kappa2=%.17g lambda_valid=%d lambda=%.17g",
						start_x, start_y, start_z, length, width, height, (int)samples.size(), (int)weights.size(), output_ess[(size_t)block], value, kappa2, lambda_valid ? 1 : 0, output_lambda[(size_t)block]);
				}
				output_values[(size_t)block] = value;
				block++;
				if(block % 5000 == 0)
				{
					append_hosc_window_tracef(
						"calculateMomentWindow.progress block=%d kept=%d last_start=%d,%d,%d last_value=%.17g",
						block,
						kept_count,
						start_x,
						start_y,
						start_z,
						value);
				}
			}
		}
	}
	append_hosc_window_tracef(
		"calculateMomentWindow.loop_end blocks=%d kept=%d expected=%d",
		block,
		kept_count,
		gridSize);
	if(writeAuditSidecar)
	{
		if(auditRunDirectory.empty() ||
		   !write_window_audit_sidecar(
			   auditRunDirectory,
			   auto_prop_name_,
			   source_grid_name_,
			   source_property_name_,
			   G_CalculateTarget,
			   al,
			   useWeightedMoment,
			   weightedOptions,
			   centerCumulantInputForOrderGE2,
			   outputStandardizedCumulant,
			   outputDiagnostics,
			   invertWindowOutput,
			   enablePercentileClipping,
			   percentileClipMin,
			   percentileClipMax,
			   varianceEpsilon,
			   runId_,
			   processId_,
			   runSequence_,
			   momentOrder,
			   momentWindowL,
			   momentWindowW,
			   momentWindowH,
			   stepX,
			   stepY,
			   is_2d ? 0 : stepZ,
			   source_origin_x_,
			   source_origin_y_,
			   source_origin_z_,
			   source_spacing_x_,
			   source_spacing_y_,
			   source_spacing_z_,
			   out_nx,
			   out_ny,
			   out_nz,
			   output_values,
			   output_support,
			   output_ess,
			   output_numeric_valid,
			   output_display_valid,
			   output_kappa2,
			   output_lambda,
			   output_lambda_valid))
		{
			append_hosc_window_trace("calculateMomentWindow.fail audit_sidecar");
			if(!auditRunDirectory.empty())
			{
				write_window_audit_status(
					auditRunDirectory,
					runId_,
					processId_,
					runSequence_,
					"failed",
					"audit_sidecar_write_failed");
			}
			return -1;
		}
		append_hosc_window_tracef("calculateMomentWindow.audit_written directory=%s", auditRunDirectory.c_str());
	}

	std::vector<double> clipped_values = output_values;
	if(enablePercentileClipping && block > 0)
	{
		std::vector<double> valid_values;
		valid_values.reserve((size_t)block);
		for(int i = 0; i < block; i++)
		{
			if(i < (int)output_keep_mask_.size() && output_keep_mask_[(size_t)i] != 0)
			{
				valid_values.push_back(output_values[(size_t)i]);
			}
		}
		if(valid_values.empty())
		{
			valid_values.assign(output_values.begin(), output_values.begin() + block);
		}
		const double lower_value = percentile_value(valid_values, percentileClipMin);
		const double upper_value = percentile_value(valid_values, percentileClipMax);
		const double clip_lower = std::min(lower_value, upper_value);
		const double clip_upper = std::max(lower_value, upper_value);
		for(int i = 0; i < block; i++)
		{
			clipped_values[(size_t)i] = clamp_double_value(output_values[(size_t)i], clip_lower, clip_upper);
		}
		append_hosc_window_tracef(
			"calculateMomentWindow.percentile_clip percentiles=%.17g,%.17g values=%.17g,%.17g",
			percentileClipMin,
			percentileClipMax,
			clip_lower,
			clip_upper);
	}

	if(output_effective_dim <= 3)
	{
		append_hosc_window_tracef(
			"calculateMomentWindow.before_create_cartesian gridSize=%d output=%d,%d,%d output_dim=%d",
			gridSize,
			out_nx,
			out_ny,
			out_nz,
			output_effective_dim);
		auto_grid_ = CreateMyGrid(&myWindowShape, "myWM");
		if(auto_grid_ == NULL)
		{
			append_hosc_window_trace("calculateMomentWindow.fail create_cartesian_null");
			return -1;
		}
		if(auto_property_ == NULL)
		{
			append_hosc_window_trace("calculateMomentWindow.fail cartesian_property_null");
			return -1;
		}
		const int output_size = auto_grid_->size();
		const int value_count = (output_size < gridSize) ? output_size : gridSize;
		for(int i = 0; i < value_count; i++)
		{
			if(i < (int)output_display_valid.size() && output_display_valid[(size_t)i] != 0)
				auto_property_->set_value(output_values[(size_t)i], i);
			else
				auto_property_->set_not_informed(i);
		}
		if(outputDiagnostics)
		{
			const string support_name = auto_prop_name_ + "_support";
			const string ess_name = auto_prop_name_ + "_ess";
			const string valid_name = auto_prop_name_ + "_numeric_valid";
			const string display_valid_name = auto_prop_name_ + "_display_valid";
			GsTLGridProperty* support_property = geostat_utils::add_property_to_grid(auto_grid_, support_name);
			GsTLGridProperty* ess_property = geostat_utils::add_property_to_grid(auto_grid_, ess_name);
			GsTLGridProperty* valid_property = geostat_utils::add_property_to_grid(auto_grid_, valid_name);
			GsTLGridProperty* display_valid_property = geostat_utils::add_property_to_grid(auto_grid_, display_valid_name);
			if(support_property == NULL || ess_property == NULL || valid_property == NULL || display_valid_property == NULL)
			{
				append_hosc_window_trace("calculateMomentWindow.fail diagnostic_property_creation");
				return -1;
			}
			for(int i = 0; i < value_count; i++)
			{
				support_property->set_value(output_support[(size_t)i], i);
				ess_property->set_value(output_ess[(size_t)i], i);
				valid_property->set_value(output_numeric_valid[(size_t)i], i);
				display_valid_property->set_value(output_display_valid[(size_t)i], i);
			}
		}
		if(outputStandardizedCumulant && G_CalculateTarget == 4)
		{
			const string kappa2_name = auto_prop_name_ + "_kappa2";
			const string lambda_name = auto_prop_name_ + "_lambda";
			const string lambda_valid_name = auto_prop_name_ + "_lambda_valid";
			GsTLGridProperty* kappa2_property = geostat_utils::add_property_to_grid(auto_grid_, kappa2_name);
			GsTLGridProperty* lambda_property = geostat_utils::add_property_to_grid(auto_grid_, lambda_name);
			GsTLGridProperty* lambda_valid_property = geostat_utils::add_property_to_grid(auto_grid_, lambda_valid_name);
			if(kappa2_property == NULL || lambda_property == NULL || lambda_valid_property == NULL)
			{
				append_hosc_window_trace("calculateMomentWindow.fail standardized_property_creation");
				return -1;
			}
			for(int i = 0; i < value_count; i++)
			{
				if(output_keep_mask_[(size_t)i] != 0 && is_float_representable(output_kappa2[(size_t)i]))
					kappa2_property->set_value(output_kappa2[(size_t)i], i);
				else
					kappa2_property->set_not_informed(i);
				if(output_lambda_valid[(size_t)i] != 0 && is_float_representable(output_lambda[(size_t)i]))
					lambda_property->set_value(output_lambda[(size_t)i], i);
				else
					lambda_property->set_not_informed(i);
				lambda_valid_property->set_value(output_lambda_valid[(size_t)i], i);
			}
		}
		if(enablePercentileClipping)
		{
			const string clipped_name = auto_prop_name_ + "_clipped";
			GsTLGridProperty* clipped_property = geostat_utils::add_property_to_grid(auto_grid_, clipped_name);
			if(clipped_property != NULL)
			{
				for(int i = 0; i < value_count; i++)
				{
					if(i < (int)output_keep_mask_.size() && output_keep_mask_[(size_t)i] != 0 &&
					   is_float_representable(clipped_values[(size_t)i]))
						clipped_property->set_value(clipped_values[(size_t)i], i);
					else
					{
						if(i < (int)output_keep_mask_.size() && output_keep_mask_[(size_t)i] != 0)
							append_hosc_window_tracef("calculateMomentWindow.clipped_float_overflow block=%d value=%.17g", i, clipped_values[(size_t)i]);
						clipped_property->set_not_informed(i);
					}
				}
			}
		}
		append_hosc_window_tracef(
			"calculateMomentWindow.after_set_cartesian_values gridSize=%d object_size=%d values=%d property=%p",
			gridSize,
			output_size,
			value_count,
			auto_property_);
		append_hosc_window_tracef("calculateMomentWindow.end blocks=%d cartesian_output=%d", block, value_count);
		if(writeAuditSidecar &&
		   !write_window_audit_status(
			   auditRunDirectory,
			   runId_,
			   processId_,
			   runSequence_,
			   "completed",
			   "host_properties_written"))
		{
			append_hosc_window_instance_trace(this, runSequence_, "audit.complete_status_failed");
			return -1;
		}
		return 0;
	}

	append_hosc_window_trace("calculateMomentWindow.fail unsupported_output_dimension");
	return -1;
}


bool HOSC_WINDOW::get_Sliding_Window(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	string wx = parameter_value_compat(parameters, "windowx.value");
	string wy = parameter_value_compat(parameters, "windowy.value");
	string wz = parameter_value_compat(parameters, "windowz.value");
	string wl = parameter_value_compat(parameters, "windowl.value");
	string ww = parameter_value_compat(parameters, "windoww.value");
	string wh = parameter_value_compat(parameters, "windowh.value");

	int iwx = atoi(wx.c_str());
	int iwy = atoi(wy.c_str());
	int iwz = atoi(wz.c_str());
	int iwl = atoi(wl.c_str());
	int iww = atoi(ww.c_str());
	int iwh = atoi(wh.c_str());

	if(iwx == -1 && iwy == -1 && iwz == -1)
		return true;
	if(iwx == 0 || iwy == 0 && iwz == 0)
	{
		error_mesgs->report("windowx","Coordinate start from 1");
	}

	if(Key->DIM == 2 && (iwz > 1 || iwh > 1))
	{
		if(iwz != 1)
		{
			error_mesgs->report("windowz","should be 1 if souce grid is 2D");
		}
		else
		{
			error_mesgs->report("windowh","should be 0/1 if souce grid is 2D");
		}
		return false;
	}



	if(iwh == 0)
		iwh = 1;
	if(iwz == 0)
		iwz = 1;


	int start = (iwz - 1) * Key->NX *Key->NY + (iwy - 1) * Key->NX + iwx;
	int totalNode = iwh * iww * iwl;
	if(totalNode == 0 || totalNode >= grid[0][0])
	{
		if(totalNode == 0)

		if(totalNode > grid[0][0])

		return true;
	}

	if(iwx + iwl - 1 > Key->NX || iwy + iww - 1 > Key->NY || (Key->DIM == 3 && iwz + iwh - 1 > Key->NZ))
	{

		return false;
	}


	if(start > grid[0][0])
	{
		error_mesgs->report("windowz","window's size out of range");
		return false;
	}

	grid[0][0] = totalNode;



	double** tempGrid = (double**)calloc(1+totalNode, sizeof(double*));
	tempGrid[0] = (double*)calloc(2, sizeof(double));
	tempGrid[0][0] = grid[0][0];
	tempGrid[0][1] = grid[0][1];

	int currNode = start;
	int currPos = 1;


	for(int w = 1; w <= totalNode; w++)
	{
		tempGrid[w] = (double*)calloc(4, sizeof(double));
		if(tempGrid[w] == NULL)
		{

			return false;
		}
		if(currNode >= 1 + (int)grid[0][0])
		{

			return false;
		}
		memcpy(tempGrid[w], grid[currNode], 4*sizeof(double));





		currNode++;
		currPos++;
		if(currPos == iwl+1)
		{
			currNode = currNode + Key->NX - iwl;
			currPos = 1;
		}
	}





	Free2DVec(&grid, 1 + (int)grid[0][0]);
	grid = tempGrid;

	Key->NX = iwl;
	Key->NY = iww;
	Key->NZ = iwh;
	Key->DIMX = iwl - 1;
	Key->DIMY = iww - 1;
	Key->DIMZ = iwh - 1;








	if(Key->NZ == 0)
	{
		Key->DIM = 2;
	}






	return true;
}




void HOSC_WINDOW::TestMyGrid()
{
	GridShape myGridShape = {10, 10, 10, 0, 0, 0, 1, 1, 1};
	RGrid* myGrid = CreateMyGrid(&myGridShape, "MyGrid");
}




RGrid* HOSC_WINDOW::CreateMyGrid(GridShape* tempGridShape, const char* gridName)
{
    GsTLlog << "HOSC_WINDOW CreateMyGrid start" << gstlIO::end;







    uintptr_t proj_addr = reinterpret_cast<uintptr_t>(current_project_);


    if (proj_addr < 0x1000) {


        return NULL;
    }


    std::string base_name(gridName);
    std::string unique_name = base_name;
    int counter = 0;

    std::string full_path;
    while (counter < 1000) {
        full_path = "/GridObject/Model/" + unique_name;
        SmartPtr<Named_interface> existing = Root::instance()->interface(full_path);

        if (!existing.raw_ptr()) {
            break;
        }

        counter++;
        std::ostringstream oss;
        oss << base_name << "_" << counter;
        unique_name = oss.str();
    }

    if (counter >= 1000) {

        return NULL;
    }


    GsTLlog << "HOSC_WINDOW CreateMyGrid before new_interface: " << full_path << gstlIO::end;
    SmartPtr<Named_interface> ni =
        Root::instance()->new_interface("cgrid", full_path);

    if (ni.raw_ptr() == 0) {

        return NULL;
    }
    GsTLlog << "HOSC_WINDOW CreateMyGrid new_interface ok" << gstlIO::end;

    Cartesian_grid* myGrid = dynamic_cast<Cartesian_grid*>(ni.raw_ptr());

    if (!myGrid) {

        return NULL;
    }
    GsTLlog << "HOSC_WINDOW CreateMyGrid cast ok" << gstlIO::end;

    int safe_nx = (tempGridShape->nx > 0) ? tempGridShape->nx : 1;
    int safe_ny = (tempGridShape->ny > 0) ? tempGridShape->ny : 1;
    int safe_nz = (tempGridShape->nz > 0) ? tempGridShape->nz : 1;
    double safe_dx = (tempGridShape->dx > 0) ? tempGridShape->dx : 1.0;
    double safe_dy = (tempGridShape->dy > 0) ? tempGridShape->dy : 1.0;
    double safe_dz = (tempGridShape->dz > 0) ? tempGridShape->dz : 1.0;
    GsTLlog << "HOSC_WINDOW CreateMyGrid dimensions request: nx=" << safe_nx << ", ny=" << safe_ny << ", nz=" << safe_nz << ", dx=" << safe_dx << ", dy=" << safe_dy << ", dz=" << safe_dz << gstlIO::end;

    myGrid->set_dimensions(safe_nx, safe_ny, safe_nz,
                          safe_dx, safe_dy, safe_dz);
    GsTLlog << "HOSC_WINDOW CreateMyGrid set_dimensions ok" << gstlIO::end;

    myGrid->origin(GsTLPoint(tempGridShape->cx, tempGridShape->cy, tempGridShape->cz));
    GsTLlog << "HOSC_WINDOW CreateMyGrid origin ok" << gstlIO::end;
    append_hosc_window_trace("CreateMyGrid: dimensions configured");



    GsTLlog << "Registering object: " << unique_name << gstlIO::end;
    append_hosc_window_tracef("CreateMyGrid.register.begin name=%s length=%d", unique_name.c_str(), (int)unique_name.size());
    unsigned long register_exception = 0;
    if(!safe_project_new_object(current_project_, &unique_name, &register_exception))
    {
        append_hosc_window_tracef("CreateMyGrid.fail safe_project_new_object exception=0x%08lx", register_exception);
        return NULL;
    }
    GsTLlog << "Grid registered successfully" << gstlIO::end;
    append_hosc_window_trace("CreateMyGrid.register.end");


    GsTLlog << "HOSC_WINDOW CreateMyGrid before managed interface lookup" << gstlIO::end;
    SmartPtr<Named_interface> managed_ni = Root::instance()->interface(full_path);
    if (!managed_ni.raw_ptr()) {
        append_hosc_window_trace("CreateMyGrid: failed to retrieve managed grid");

        return NULL;
    }
    GsTLlog << "HOSC_WINDOW CreateMyGrid managed interface ok" << gstlIO::end;

    RGrid* managed_grid = dynamic_cast<RGrid*>(managed_ni.raw_ptr());
    Cartesian_grid* managed_cgrid = dynamic_cast<Cartesian_grid*>(managed_ni.raw_ptr());
    if(managed_grid == NULL || managed_cgrid == NULL)
    {
        append_hosc_window_trace("CreateMyGrid: managed grid cast failed");
        return NULL;
    }

    GsTLlog << "HOSC_WINDOW CreateMyGrid before add_property_to_grid: " << auto_prop_name_ << gstlIO::end;
	auto_property_ = geostat_utils::add_property_to_grid(managed_grid, auto_prop_name_);
    if(auto_property_ == NULL)
    {
        append_hosc_window_trace("CreateMyGrid: add_property_to_grid failed");
        return NULL;
    }
    GsTLlog << "HOSC_WINDOW CreateMyGrid add_property_to_grid ok" << gstlIO::end;
	managed_grid->set_level(1);
	auto_grid_name_ = unique_name;
	auto_grid_registered_ = true;
    append_hosc_window_trace("CreateMyGrid: property attached to managed grid");

    return managed_grid;
}

Geostat_grid* HOSC_WINDOW::CreateSparseOutputPointSet(
	const char* gridName,
	const std::vector<Geostat_grid::location_type>& point_locations,
	const std::vector<double>& point_values)
{
	GsTLlog << "HOSC_WINDOW CreateSparseOutputPointSet start" << gstlIO::end;
	append_hosc_window_tracef(
		"CreateSparseOutputPointSet.begin points=%d values=%d",
		(int)point_locations.size(),
		(int)point_values.size());
	uintptr_t proj_addr = reinterpret_cast<uintptr_t>(current_project_);
	if(proj_addr < 0x1000)
	{
		append_hosc_window_trace("CreateSparseOutputPointSet.fail invalid_project");
		return NULL;
	}
	if(point_locations.size() != point_values.size())
	{
		append_hosc_window_trace("CreateSparseOutputPointSet.fail size_mismatch");
		return NULL;
	}

	std::string base_name(gridName);
	std::string unique_name = base_name;
	int counter = 0;
	std::string full_path;
	while(counter < 1000)
	{
		full_path = "/GridObject/Model/" + unique_name;
		SmartPtr<Named_interface> existing = Root::instance()->interface(full_path);
		if(!existing.raw_ptr())
		{
			break;
		}
		counter++;
		std::ostringstream oss;
		oss << base_name << "_" << counter;
		unique_name = oss.str();
	}
	if(counter >= 1000)
	{
		append_hosc_window_trace("CreateSparseOutputPointSet.fail unique_name");
		return NULL;
	}

	char temp_directory[MAX_PATH] = "";
	char temp_file[MAX_PATH] = "";
	const DWORD temp_length = GetTempPathA(MAX_PATH, temp_directory);
	if(temp_length == 0 || temp_length >= MAX_PATH ||
	   GetTempFileNameA(temp_directory, "hsc", 0, temp_file) == 0)
	{
		append_hosc_window_trace("CreateSparseOutputPointSet.fail temporary_path");
		return NULL;
	}
	const std::string tmp_filename(temp_file);
	if(!write_sgems_pointset_file(tmp_filename, unique_name, auto_prop_name_, point_locations, point_values))
	{
		append_hosc_window_trace("CreateSparseOutputPointSet.fail write_sgems");
		return NULL;
	}

	Error_messages_handler load_errors;
	const std::string params = tmp_filename + "::s-gems";
	append_hosc_window_tracef(
		"CreateSparseOutputPointSet.before_load action=LoadObjectFromFile params=%s",
		params.c_str());
	unsigned long load_exception = 0;
	const std::string load_action = "LoadObjectFromFile";
	const bool load_ok = safe_project_execute_action(
		current_project_,
		&load_action,
		&params,
		&load_errors,
		&load_exception);
	if(!load_ok)
	{
		append_hosc_window_tracef(
			"CreateSparseOutputPointSet.warn load_object_failed exception=0x%08lx errors=%s",
			load_exception,
			load_errors.errors().c_str());
	}

	SmartPtr<Named_interface> loaded_ni = Root::instance()->interface(full_path);
	Geostat_grid* loaded_grid = dynamic_cast<Geostat_grid*>(loaded_ni.raw_ptr());
	if(loaded_grid == NULL)
	{
		append_hosc_window_tracef("CreateSparseOutputPointSet.fail loaded_grid_null path=%s", full_path.c_str());
		return NULL;
	}
	auto_property_ = loaded_grid->property(auto_prop_name_);
	if(auto_property_ == NULL)
	{
		append_hosc_window_trace("CreateSparseOutputPointSet.fail loaded_property_null");
		return NULL;
	}

	auto_grid_registered_ = true;
	auto_grid_name_ = unique_name;
	append_hosc_window_tracef(
		"CreateSparseOutputPointSet.end name=%s points=%d property=%s",
		unique_name.c_str(),
		(int)point_locations.size(),
		auto_prop_name_.c_str());
	return loaded_grid;
}

int HOSC_WINDOW::set_MyGrid_Data(double** targetGrid, int gridSize, int valuePos)
{


	for(int i = 0; i < gridSize; i++)
	{
		auto_property_->set_value(targetGrid[i][valuePos], i);

	}


	return 0;
}

bool HOSC_WINDOW::get_My_Shape(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	auto_prop_name_ = parameter_value_compat(parameters, "props.value" );
	if(auto_prop_name_.empty())
	{
		error_mesgs->report("props","No property specified");
		return false;
	}


	string scx = parameter_value_compat(parameters, "xmn_.value" );
	string scy = parameter_value_compat(parameters, "ymn_.value" );
	string scz = parameter_value_compat(parameters, "zmn_.value" );





	double dcx = atof(scx.c_str());
	double dcy = atof(scy.c_str());
	double dcz = atof(scz.c_str());




	string snx = parameter_value_compat(parameters, "fst_dir.value" );
	string sny = parameter_value_compat(parameters, "snd_dir.value" );
	string snz = parameter_value_compat(parameters, "trd_dir.value" );

	int inx = atoi(snx.c_str());
	int iny = atoi(sny.c_str());
	int inz = atoi(snz.c_str());



	if(inx > Key->NDIR || iny > Key->NDIR || inz > Key->NDIR || inz < 0 || iny < 0 || inx < 0 || inx == iny || inx == inz || iny == inz)
	{
		error_mesgs->report("fst_dir","Out of range");
		error_mesgs->report("snd_dir","Out of range");
		error_mesgs->report("trd_dir","Out of range");
		return false;
	}

	int iinx = (inx == 0) ? 1:Key->NLAG[inx-1];
	int iiny = (iny == 0) ? 1:Key->NLAG[iny-1];
	int iinz = (inz == 0) ? 1:Key->NLAG[inz-1];

	int dirNum = 3;
	if(iinx == 1 || iiny == 1 || iinz == 1)
		dirNum = 2;

	GridShape tempShape = {iinx, iiny, iinz, dcx, dcy, dcz, 1, 1, 1};

	string crossectionContent = parameter_value_compat(parameters, "cross_pos.value" );
	char temp[256];
	strcpy(temp, crossectionContent.c_str());
	int count = 0;
    int len = strlen(temp);

    if (len == 0) {
        error_mesgs->report("cross_pos","Wrong input");
        return false;
    }

    if (temp[0] == ',' || temp[len-1] == ',') {
        error_mesgs->report("cross_pos","Wrong input");
        return false;
    }

    for (int i = 0; i < len; i++) {
        if (!(isdigit(temp[i]) || temp[i] == ',' || temp[i] == '-')) {
            error_mesgs->report("cross_pos","Wrong input");
            return false;
        }
        if (temp[i] == ',' && temp[i+1] == ',') {
            error_mesgs->report("cross_pos","Wrong input");
            return false;
        }
    }
	crossectionNum = 0;
	char *token = strtok(temp, ",");
    while (token != NULL) {
        crossection[crossectionNum++] = atoi(token);
		if(crossection[crossectionNum-1] <= 0)
		{
			error_mesgs->report("cross_pos","Wrong input");
			return false;
		}
        token = strtok(NULL, ",");
    }

	if(Key->NDIR - dirNum == 0)
	{
		if(crossectionNum > 1)
		{
			error_mesgs->report("cross_pos", "too many coordinates for this crossection");
			return false;
		}
		crossection[0] = 1;
	}
	else if(crossectionNum != Key->NDIR - dirNum)
	{
		error_mesgs->report("cross_pos","coordinates num should equal to Key->NDIR - dirNum");
		return false;
	}







	myGridShape = tempShape;



	int* tempNOD = (int*)calloc(Key->NDIR, sizeof(int));
	int* tempNLAG = (int*)calloc(Key->NDIR, sizeof(int));
	int* tempLAG = (int*)calloc(Key->NDIR, sizeof(int));
	int* tempAZM = (int*)calloc(Key->NDIR, sizeof(int));
	int* tempATOL = (int*)calloc(Key->NDIR, sizeof(int));
	int* tempBANDWH = (int*)calloc(Key->NDIR, sizeof(int));
	int* tempDIP = (int*)calloc(Key->NDIR, sizeof(int));
	int* tempDTOL = (int*)calloc(Key->NDIR, sizeof(int));
	int* tempBANDWD = (int*)calloc(Key->NDIR, sizeof(int));

	if(tempNLAG == NULL || tempLAG == NULL || tempAZM == NULL || tempATOL == NULL || tempBANDWH == NULL || tempDIP == NULL || tempDTOL == NULL || tempBANDWD == NULL || tempNOD == NULL)
	{

		return false;
	}
	for(int i = 0; i < Key->NDIR - 3; i++)
	{
		tempNLAG[i] = 1;
	}

	if(iinx == 1)
	{
		int tempI = 0;
		for(int i = 0; i < Key->NDIR; i++)
		{
			if(i==iny-1 || i==inz-1)
				continue;
			tempNOD[tempI] = Key->NOD[i];
			tempNLAG[tempI] = Key->NLAG[i];
			tempLAG[tempI] = Key->LAG[i];
			tempAZM[tempI] = Key->AZM[i];
			tempATOL[tempI] = Key->ATOL[i];
			tempBANDWH[tempI] = Key->BANDWD[i];
			tempDIP[tempI] = Key->DIP[i];
			tempDTOL[tempI] = Key->DTOL[i];
			tempBANDWD[tempI] = Key->BANDWD[i];
			tempI++;
		}
		tempNLAG[Key->NDIR-1] = iiny;
		tempNLAG[Key->NDIR-2] = iinz;

		tempNOD[Key->NDIR-1] = Key->NOD[iny-1];
		tempNOD[Key->NDIR-2] = Key->NOD[inz-1];
		tempLAG[Key->NDIR-1] = Key->LAG[iny-1];
		tempLAG[Key->NDIR-2] = Key->LAG[inz-1];
		tempAZM[Key->NDIR-1] = Key->AZM[iny-1];
		tempAZM[Key->NDIR-2] = Key->AZM[inz-1];
		tempATOL[Key->NDIR-1] = Key->ATOL[iny-1];
		tempATOL[Key->NDIR-2] = Key->ATOL[inz-1];
		tempBANDWH[Key->NDIR-1] = Key->BANDWH[iny-1];
		tempBANDWH[Key->NDIR-2] = Key->BANDWH[inz-1];
		tempDIP[Key->NDIR-1] = Key->DIP[iny-1];
		tempDIP[Key->NDIR-2] = Key->DIP[inz-1];
		tempDTOL[Key->NDIR-1] = Key->DTOL[iny-1];
		tempDTOL[Key->NDIR-2] = Key->DTOL[inz-1];
		tempBANDWD[Key->NDIR-1] = Key->BANDWD[iny-1];
		tempBANDWD[Key->NDIR-2] = Key->BANDWD[inz-1];
	}
	else if(iiny == 1)
	{
		int tempI = 0;
		for(int i = 0; i < Key->NDIR; i++)
		{
			if(i==inx-1 || i==inz-1)
				continue;
			tempNOD[tempI] = Key->NOD[i];
			tempNLAG[tempI] = Key->NLAG[i];
			tempLAG[tempI] = Key->LAG[i];
			tempAZM[tempI] = Key->AZM[i];
			tempATOL[tempI] = Key->ATOL[i];
			tempBANDWH[tempI] = Key->BANDWD[i];
			tempDIP[tempI] = Key->DIP[i];
			tempDTOL[tempI] = Key->DTOL[i];
			tempBANDWD[tempI] = Key->BANDWD[i];
			tempI++;
		}
		tempNLAG[Key->NDIR-1] = iinx;
		tempNLAG[Key->NDIR-2] = iinz;

		tempNOD[Key->NDIR-1] = Key->NOD[inx-1];
		tempNOD[Key->NDIR-2] = Key->NOD[inz-1];
		tempLAG[Key->NDIR-1] = Key->LAG[inx-1];
		tempLAG[Key->NDIR-2] = Key->LAG[inz-1];
		tempAZM[Key->NDIR-1] = Key->AZM[inx-1];
		tempAZM[Key->NDIR-2] = Key->AZM[inz-1];
		tempATOL[Key->NDIR-1] = Key->ATOL[inx-1];
		tempATOL[Key->NDIR-2] = Key->ATOL[inz-1];
		tempBANDWH[Key->NDIR-1] = Key->BANDWH[inx-1];
		tempBANDWH[Key->NDIR-2] = Key->BANDWH[inz-1];
		tempDIP[Key->NDIR-1] = Key->DIP[inx-1];
		tempDIP[Key->NDIR-2] = Key->DIP[inz-1];
		tempDTOL[Key->NDIR-1] = Key->DTOL[inx-1];
		tempDTOL[Key->NDIR-2] = Key->DTOL[inz-1];
		tempBANDWD[Key->NDIR-1] = Key->BANDWD[inx-1];
		tempBANDWD[Key->NDIR-2] = Key->BANDWD[inz-1];
	}
	else if(iinz == 1)
	{
		int tempI = 0;
		for(int i = 0; i < Key->NDIR; i++)
		{
			if(i==inx-1 || i==iny-1)
				continue;
			tempNOD[tempI] = Key->NOD[i];
			tempNLAG[tempI] = Key->NLAG[i];
			tempLAG[tempI] = Key->LAG[i];
			tempAZM[tempI] = Key->AZM[i];
			tempATOL[tempI] = Key->ATOL[i];
			tempBANDWH[tempI] = Key->BANDWD[i];
			tempDIP[tempI] = Key->DIP[i];
			tempDTOL[tempI] = Key->DTOL[i];
			tempBANDWD[tempI] = Key->BANDWD[i];
			tempI++;
		}
		tempNLAG[Key->NDIR-1] = iinx;
		tempNLAG[Key->NDIR-2] = iiny;

		tempNOD[Key->NDIR-1] = Key->NOD[inx-1];
		tempNOD[Key->NDIR-2] = Key->NOD[iny-1];
		tempLAG[Key->NDIR-1] = Key->LAG[inx-1];
		tempLAG[Key->NDIR-2] = Key->LAG[iny-1];
		tempAZM[Key->NDIR-1] = Key->AZM[inx-1];
		tempAZM[Key->NDIR-2] = Key->AZM[iny-1];
		tempATOL[Key->NDIR-1] = Key->ATOL[inx-1];
		tempATOL[Key->NDIR-2] = Key->ATOL[iny-1];
		tempBANDWH[Key->NDIR-1] = Key->BANDWH[inx-1];
		tempBANDWH[Key->NDIR-2] = Key->BANDWH[iny-1];
		tempDIP[Key->NDIR-1] = Key->DIP[inx-1];
		tempDIP[Key->NDIR-2] = Key->DIP[iny-1];
		tempDTOL[Key->NDIR-1] = Key->DTOL[inx-1];
		tempDTOL[Key->NDIR-2] = Key->DTOL[iny-1];
		tempBANDWD[Key->NDIR-1] = Key->BANDWD[inx-1];
		tempBANDWD[Key->NDIR-2] = Key->BANDWD[iny-1];
	}
	else if(iinx != 1 && iiny != 1 && iinz != 1)
	{
		int tempI = 0;
		for(int i = 0; i < Key->NDIR; i++)
		{
			if(i==inx-1 ||i==iny-1 || i==inz-1)
				continue;
			tempNOD[tempI] = Key->NOD[i];
			tempNLAG[tempI] = Key->NLAG[i];
			tempLAG[tempI] = Key->LAG[i];
			tempAZM[tempI] = Key->AZM[i];
			tempATOL[tempI] = Key->ATOL[i];
			tempBANDWH[tempI] = Key->BANDWD[i];
			tempDIP[tempI] = Key->DIP[i];
			tempDTOL[tempI] = Key->DTOL[i];
			tempBANDWD[tempI] = Key->BANDWD[i];
			tempI++;
		}
		tempNLAG[Key->NDIR-1] = iinx;
		tempNLAG[Key->NDIR-2] = iiny;
		tempNLAG[Key->NDIR-3] = iinz;
		tempNOD[Key->NDIR-1] = Key->NOD[inx-1];
		tempNOD[Key->NDIR-2] = Key->NOD[iny-1];
		tempNOD[Key->NDIR-3] = Key->NOD[inz-1];
		tempLAG[Key->NDIR-1] = Key->LAG[inx-1];
		tempLAG[Key->NDIR-2] = Key->LAG[iny-1];
		tempLAG[Key->NDIR-3] = Key->LAG[inz-1];
		tempAZM[Key->NDIR-1] = Key->AZM[inx-1];
		tempAZM[Key->NDIR-2] = Key->AZM[iny-1];
		tempAZM[Key->NDIR-3] = Key->AZM[inz-1];
		tempATOL[Key->NDIR-1] = Key->ATOL[inx-1];
		tempATOL[Key->NDIR-2] = Key->ATOL[iny-1];
		tempATOL[Key->NDIR-3] = Key->ATOL[inz-1];
		tempBANDWH[Key->NDIR-1] = Key->BANDWH[inx-1];
		tempBANDWH[Key->NDIR-2] = Key->BANDWH[iny-1];
		tempBANDWH[Key->NDIR-3] = Key->BANDWH[inz-1];
		tempDIP[Key->NDIR-1] = Key->DIP[inx-1];
		tempDIP[Key->NDIR-2] = Key->DIP[iny-1];
		tempDIP[Key->NDIR-3] = Key->DIP[inz-1];
		tempDTOL[Key->NDIR-1] = Key->DTOL[inx-1];
		tempDTOL[Key->NDIR-2] = Key->DTOL[iny-1];
		tempDTOL[Key->NDIR-3] = Key->DTOL[inz-1];
		tempBANDWD[Key->NDIR-1] = Key->BANDWD[inx-1];
		tempBANDWD[Key->NDIR-2] = Key->BANDWD[iny-1];
		tempBANDWD[Key->NDIR-3] = Key->BANDWD[inz-1];
	}
	for(int i = 0; i < Key->NDIR; i++)
	{
		Key->NOD[i] = tempNOD[i];
		Key->NLAG[i] = tempNLAG[i];
		Key->LAG[i] = tempLAG[i];
		Key->AZM[i] = tempAZM[i];
		Key->ATOL[i] = tempATOL[i];
		Key->BANDWH[i] = tempBANDWH[i];
		Key->DIP[i] = tempDIP[i];
		Key->DTOL[i] = tempDTOL[i];
		Key->BANDWD[i] = tempBANDWD[i];
	}

	if(Key->NDIR - dirNum != 0)
		for(int i = 0; i < crossectionNum; i++)
		{
			if(crossection[i] > Key->NLAG[i])
			{
				error_mesgs->report("cross_pos","coordinates out of range");
				return false;
			}
		}

	return true;

}

#if 0
RGrid* HOSC_WINDOW::getDefaultGrid()
{
    std::string grid_path = gridModels_manager + "/HOSC_DefaultGrid";
    SmartPtr<Named_interface> ni = Root::instance()->interface(grid_path);
    return dynamic_cast<RGrid*>(ni.raw_ptr());
}

bool HOSC_WINDOW::get_Default_Data(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	if(default_grid_ == NULL)
		return false;
	string default_prop_name_ = parameter_value_compat(parameters, "props.value" );
	if(default_prop_name_.empty())
	{
		error_mesgs->report("props","No property specified");
		return false;
	}
	default_property_name_ = geostat_utils::add_property_to_grid(default_grid_, default_prop_name_);

	GsTLGridProperty** KeyProp = (GsTLGridProperty**)calloc(Key->NDIR, sizeof(GsTLGridProperty*));
	if(KeyProp == NULL)
	{

		return NULL;
	}

    for(int i = 0; i < Key->NDIR; i++)
	{
		char c = 'h';
		char result[32];
		sprintf(result, "%c%d", c, i);
		KeyProp[i] = geostat_utils::add_property_to_grid(default_grid_, result);
	}
	return true;
}
#endif

#if 0
namespace {
    class HOSC_Init {
    public:
        HOSC_Init() {
            GsTLlog << "HOSC_WINDOW Plugin initializing..." << gstlIO::end;


            SmartPtr<Named_interface> ui_ni =
                Root::instance()->new_interface(
                    "HOSC_WINDOW.ui",
                    "HOSC_WINDOW"
                );

            if (ui_ni.raw_ptr()) {
                GsTLlog << "HOSC_WINDOW UI registered successfully" << gstlIO::end;
            }


            createDefaultGrid();
        }

    private:
        void createDefaultGrid() {

            std::string full_path = "/GridObject/Model/HOSC_DefaultGrid";


            SmartPtr<Named_interface> existing = Root::instance()->interface(full_path);
            if (existing.raw_ptr()) {
                GsTLlog << "Default grid already exists" << gstlIO::end;
                return;
            }


            SmartPtr<Named_interface> ni =
                Root::instance()->new_interface("cgrid", full_path);

            if (!ni.raw_ptr()) {
                GsTLlog << "Failed to create default grid" << gstlIO::end;
                return;
            }

            Cartesian_grid* grid = dynamic_cast<Cartesian_grid*>(ni.raw_ptr());
            if (grid) {

                grid->set_dimensions(0, 0, 0, 0, 0, 0);
                grid->origin(GsTLPoint(0, 0, 0));


                SmartPtr<Named_interface> proj_ni =
                    Root::instance()->interface(projects_manager);

                if (proj_ni.raw_ptr()) {
                    Manager* proj_mgr = dynamic_cast<Manager*>(proj_ni.raw_ptr());
                    if (proj_mgr) {
                        Manager::interface_iterator it = proj_mgr->begin_interfaces();
                        if (it != proj_mgr->end_interfaces()) {
                            GsTL_project* project = dynamic_cast<GsTL_project*>((*it).raw_ptr());
                            if (project) {
                                project->new_object("HOSC_DefaultGrid");
                                GsTLlog << "Default grid created and registered" << gstlIO::end;
                            }
                        }
                    }
                }
            }
        }
    };

    static HOSC_Init hosc_initializer;
}

#endif


GEOSTAT_PLUGIN(HOSC_WINDOW);
