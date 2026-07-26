
#define NOMINMAX
#include <algorithm>
#include <cmath>
#include <GsTLAppli/geostat/parameters_handler.h>
#include <GsTLAppli/geostat/parameters_handler_impl.h>
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
#include <limits.h>
#include <new>
#include <stdarg.h>
#include <sstream>
#include <fstream>
#include <time.h>
#include <excpt.h>

#include "Plugin.h"
#include "PCA.h"

namespace wmc = weighted_moment_core;

namespace {
std::string parameter_value_compat(const Parameters_handler* parameters, const char* name)
{
    if(parameters == NULL || name == NULL)
    {
        return std::string();
    }






    struct Parameters_handler_xml_layout
    {
        void* vptr;
        long references;
        bool is_ready;
        unsigned char padding[3];
        QDomDocument document;
        QDomElement root;
    };

    const Parameters_handler_xml_layout* xml =
        reinterpret_cast<const Parameters_handler_xml_layout*>(parameters);
    const char* dot = strchr(name, '.');
    if(dot == NULL || dot == name || dot[1] == '\0')
    {
        return std::string();
    }

    char element_name[128];
    const size_t element_length = (size_t)(dot - name);
    if(element_length >= sizeof(element_name))
    {
        return std::string();
    }
    memcpy(element_name, name, element_length);
    element_name[element_length] = '\0';

    const QDomNodeList nodes =
        xml->document.elementsByTagName(QString(element_name));
    if(nodes.count() <= 0 || !nodes.item(0).isElement())
    {
        return std::string();
    }
    const QDomElement element = nodes.item(0).toElement();
    const QString value = element.attribute(QString(dot + 1));
    const char* latin = value.latin1();
    return latin == NULL ? std::string() : std::string(latin);
}

void append_hosc_window_trace(const char* message)
{
    FILE* trace = fopen("HOSC_WINDOW_debug.log", "a");
    if(trace == NULL)
    {
        return;
    }
    fprintf(trace, "%ld %s\n", (long)time(NULL), message == NULL ? "" : message);
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

wmc::WeightMode parse_weight_mode_text(const std::string& text)
{
    const std::string normalized = normalize_algorithm_text(text);
    if(normalized == "coordinate_gaussian" || normalized == "coordinate gaussian" || normalized == "spatial_gaussian" || normalized == "spatial gaussian" || normalized == "coordinate" || normalized == "coord" || normalized == "coord_gauss")
    {
        return wmc::WeightMode_CoordinateGaussian;
    }
    if(normalized == "value_gaussian" || normalized == "value gaussian" || normalized == "value" || normalized == "value_gauss")
    {
        return wmc::WeightMode_ValueGaussian;
    }
    return wmc::WeightMode_Uniform;
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
    case wmc::WeightMode_ValueGaussian: return "WVG";
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
    __except(EXCEPTION_EXECUTE_HANDLER)
    {
        if(exception_code != NULL)
        {
            *exception_code = 1;
        }
        append_hosc_window_trace("safe_calculateMomentWindow.exception");
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
    return grid_obj->nx() <= 1 || grid_obj->ny() <= 1 || grid_obj->nz() <= 1;
}

bool is_effective_2d_dims(int nx, int ny, int nz)
{
    return nx <= 1 || ny <= 1 || nz <= 1;
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
	append_hosc_window_trace("HOSC_WINDOW.ctor build=movewindow_full_window_cartesian_20260726");
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
	percentileClipMin = 2.0;
	percentileClipMax = 98.0;
	useWeightedMoment = false;
	weightedOptions = wmc::WeightedMomentOptions();

	G_CalculateTarget = -1;
	G_TargetDimension = 0;




}

HOSC_WINDOW::~HOSC_WINDOW()
{
	append_hosc_window_tracef(
		"HOSC_WINDOW.dtor.begin this=%p grid=%p rows=%d source_grid=%p auto_grid=%p auto_property=%p",
		this,
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
	append_hosc_window_tracef("HOSC_WINDOW.dtor.end this=%p", this);
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
		else if(weightedOptions.mode == wmc::WeightMode_ValueGaussian)
		{
			weight_token = 'V';
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
	GsTLlog << "HOSC_WINDOW initialize start" << gstlIO::end;
	append_hosc_window_trace("initialize.begin");
	get_Calculate_Target(parameters, error_mesgs);
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
	GsTLlog << "HOSC_WINDOW source ready: nx=" << source_nx_ << ", ny=" << source_ny_ << ", nz=" << source_nz_ << gstlIO::end;
	append_hosc_window_tracef(
		"initialize.source_ready nx=%d ny=%d nz=%d nxyz=%d raw_nxyz=%d grid=%p source_val=%p",
		source_nx_,
		source_ny_,
		source_nz_,
		(int)grid[0][0],
		source_grid_->size(),
		source_grid_,
		source_val);











	if( !get_Windows_Moment(parameters, error_mesgs))
	{
		append_hosc_window_trace("initialize.fail windows_moment");
		return false;
	}
	GsTLlog << "HOSC_WINDOW window ready: wx=" << momentWindowL << ", wy=" << momentWindowW << ", wz=" << momentWindowH << ", order=" << momentOrder << ", out_nx=" << myWindowShape.nx << ", out_ny=" << myWindowShape.ny << ", out_nz=" << myWindowShape.nz << gstlIO::end;
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
		 FILE* file = fopen("tempData.txt", "w");
		if(file == NULL)
		{
			append_hosc_window_trace("execute.warn tempData_open_failed");
			return 0;
		}
		fprintf(file, "Order%d\n", momentOrder);
		fprintf(file, "%d\n", 5);
		fprintf(file, "index\n");
		fprintf(file, "x\ny\nz\nvalue\n");
		for(int i = 0; i < auto_grid_->size(); i++)
		{
			fprintf(file, "%d    ", i);
			Geostat_grid::location_type loc = auto_grid_->location(i);
			fprintf(file, "%.15lf    ", (double)loc.x());
			fprintf(file, "%.15lf    ", (double)loc.y());
			fprintf(file, "%.15lf    ", (double)loc.z());
			fprintf(file, "%.15lf\n", (double)auto_property_->get_value(i));
		}
		fclose(file);
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
    const bool source_is_regular = source_rgrid_ != NULL;
    if(!source_grid_)
    {
        std::ostringstream error_stream;
        error_stream << grid_name << " is not a valid geostatistical grid object";
        error_mesgs->report("Source_Grid", error_stream.str());
        append_hosc_window_trace("get_Source_Data.fail invalid_geostat_grid");
        return false;
    }
    append_hosc_window_tracef(
        "get_Source_Data.grid class=%s is_regular=%d raw_size=%d",
        source_grid_->classname().c_str(),
        source_is_regular ? 1 : 0,
        source_grid_->size());

    source_x = source_grid_->select_property("x");
    if(source_x == NULL)
        source_x = source_grid_->select_property("X");
    source_y = source_grid_->select_property("y");
    if(source_y == NULL)
        source_y = source_grid_->select_property("Y");
    source_z = source_grid_->select_property("z");
    if(source_z == NULL)
        source_z = source_grid_->select_property("Z");

    string prop_name = parameter_value_compat(parameters, "Source_Property.value");
    append_hosc_window_tracef("get_Source_Data.props requested=%s has_x=%d has_y=%d has_z=%d", prop_name.c_str(), source_x != NULL, source_y != NULL, source_z != NULL);
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
    if(prop_name == "z" || prop_name == "Z")
    {
        source_z = NULL;
        append_hosc_window_trace("get_Source_Data.info source property z treated as value column; spatial z derived from regular grid index");
    }

    const int input_count = source_grid_->size();
    if(input_count <= 0)
    {
        error_mesgs->report("Source_Grid", "Source grid has no nodes");
        append_hosc_window_trace("get_Source_Data.fail empty_source_grid");
        return false;
    }

    const bool has_explicit_coords = source_x != NULL && source_y != NULL && source_z != NULL;
	Cartesian_grid* source_cartesian = dynamic_cast<Cartesian_grid*>(source_rgrid_);
	if(source_cartesian != NULL && !has_explicit_coords)
	{
		const GsTLPoint source_origin = source_cartesian->origin();
		const GsTLCoordVector source_spacing = source_cartesian->cell_dimensions();
		source_origin_x_ = source_origin.x();
		source_origin_y_ = source_origin.y();
		source_origin_z_ = source_origin.z();
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
    const int safe_raw_nx = (source_rgrid_ != NULL && source_rgrid_->nx() > 0) ? source_rgrid_->nx() : 1;
    const int safe_raw_ny = (source_rgrid_ != NULL && source_rgrid_->ny() > 0) ? source_rgrid_->ny() : 1;
    const int raw_xy = safe_raw_nx * safe_raw_ny;

    for(int i = 0; i < input_count; i++)
    {
        const int zero_index = i;
        int x = zero_index % safe_raw_nx + 1;
        int y = (zero_index / safe_raw_nx) % safe_raw_ny + 1;
        int z = raw_xy > 0 ? zero_index / raw_xy + 1 : 1;
        if(!source_is_regular)
        {
            const Geostat_grid::location_type loc = source_grid_->location(i);
            x = round_grid_coordinate(loc.x());
            y = round_grid_coordinate(loc.y());
            z = round_grid_coordinate(loc.z());
        }
        else if(has_explicit_coords)
        {
            x = (int)source_x->get_value(i);
            y = (int)source_y->get_value(i);
            z = (int)source_z->get_value(i);
        }
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
	if(source_cartesian == NULL || has_explicit_coords)
	{
		source_origin_x_ = (double)min_x;
		source_origin_y_ = (double)min_y;
		source_origin_z_ = (double)min_z;
		source_spacing_x_ = 1.0;
		source_spacing_y_ = 1.0;
		source_spacing_z_ = 1.0;
	}
    const int dense_count = source_nx_ * source_ny_ * source_nz_;
    if(source_nx_ <= 0 || source_ny_ <= 0 || source_nz_ <= 0 || dense_count <= 0)
    {
        error_mesgs->report("Source_Grid", "Inferred source grid dimensions are invalid");
        append_hosc_window_tracef("get_Source_Data.fail invalid_dense_dims dims=%d,%d,%d", source_nx_, source_ny_, source_nz_);
        return false;
    }

    append_hosc_window_tracef(
        "get_Source_Data.densify inferred=%d,%d,%d dense_count=%d input_count=%d min=%d,%d,%d max=%d,%d,%d explicit_coords=%d non_regular=%d",
        source_nx_, source_ny_, source_nz_, dense_count, input_count,
        min_x, min_y, min_z, max_x, max_y, max_z,
        has_explicit_coords ? 1 : 0,
        source_is_regular ? 0 : 1);

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
        int raw_x = zero_index % safe_raw_nx + 1;
        int raw_y = (zero_index / safe_raw_nx) % safe_raw_ny + 1;
        int raw_z = raw_xy > 0 ? zero_index / raw_xy + 1 : 1;
        if(!source_is_regular)
        {
            const Geostat_grid::location_type loc = source_grid_->location(i);
            raw_x = round_grid_coordinate(loc.x());
            raw_y = round_grid_coordinate(loc.y());
            raw_z = round_grid_coordinate(loc.z());
        }
        else if(has_explicit_coords)
        {
            raw_x = (int)source_x->get_value(i);
            raw_y = (int)source_y->get_value(i);
            raw_z = (int)source_z->get_value(i);
        }
        const int local_x = raw_x - min_x;
        const int local_y = raw_y - min_y;
        const int local_z = raw_z - min_z;
        if(local_x < 0 || local_x >= source_nx_ || local_y < 0 || local_y >= source_ny_ || local_z < 0 || local_z >= source_nz_)
        {
            append_hosc_window_tracef("get_Source_Data.skip_oob raw=%d,%d,%d local=%d,%d,%d", raw_x, raw_y, raw_z, local_x, local_y, local_z);
            continue;
        }
        const int dense_index = local_z * source_nx_ * source_ny_ + local_y * source_nx_ + local_x + 1;
        grid[dense_index][3] = (double)source_val->get_value(i);
        source_present_[(size_t)(dense_index - 1)] = 1;
        filled_count++;
    }
	if(filled_count == 0)
	{
		error_mesgs->report("Source_Property", "Selected source property has no informed values");
		append_hosc_window_trace("get_Source_Data.fail no_informed_values");
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
		append_hosc_window_trace("get_Calculate_Target.warn unknown_target");

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
	else
		al = AVL;
	append_hosc_window_tracef("get_Calculate_Algorithm.parsed algorithm=%d", (int)al);
	append_hosc_window_trace("get_Calculate_Algorithm.before_return");
	return true;
}

bool HOSC_WINDOW::get_Weighted_Moment_Options(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	append_hosc_window_trace("get_Weighted_Moment_Options.begin");
	weightedOptions = wmc::WeightedMomentOptions();
	useWeightedMoment = false;

	string enabled = parameter_value_compat(parameters, "enable_weighted_moment.value");
	if(enabled.empty())
	{
		enabled = parameter_value_compat(parameters, "Use_Weighted_Moment.value");
	}
	useWeightedMoment = parse_bool_text(enabled, false);

	const string mode = parameter_value_compat(parameters, "weight_mode.value");
	if(!mode.empty())
	{
		weightedOptions.mode = parse_weight_mode_text(mode);
	}

	const string bandwidth = parameter_value_compat(parameters, "weight_bandwidth.value");
	if(!bandwidth.empty())
	{
		weightedOptions.bandwidth = atof(bandwidth.c_str());
	}

	const string reference_value = parameter_value_compat(parameters, "weight_reference_value.value");
	if(!reference_value.empty())
	{
		weightedOptions.reference_value = atof(reference_value.c_str());
	}

	const string center_x = parameter_value_compat(parameters, "weight_center_x.value");
	const string center_y = parameter_value_compat(parameters, "weight_center_y.value");
	const string center_z = parameter_value_compat(parameters, "weight_center_z.value");
	if(!center_x.empty()) weightedOptions.center_x = atof(center_x.c_str());
	if(!center_y.empty()) weightedOptions.center_y = atof(center_y.c_str());
	if(!center_z.empty()) weightedOptions.center_z = atof(center_z.c_str());

	const string min_weight = parameter_value_compat(parameters, "min_weight.value");
	if(!min_weight.empty())
	{
		weightedOptions.min_weight = atof(min_weight.c_str());
	}

	const string normalize = parameter_value_compat(parameters, "normalize_weights.value");
	if(!normalize.empty())
	{
		weightedOptions.normalize_weights = parse_bool_text(normalize, true);
	}

	const string zero_first = parameter_value_compat(parameters, "zero_first_order_moment.value");
	if(!zero_first.empty())
	{
		weightedOptions.zero_first_order_moment = parse_bool_text(zero_first, false);
	}

	if(weightedOptions.bandwidth <= 0.0)
	{
		append_hosc_window_trace("get_Weighted_Moment_Options.warn bandwidth<=0 reset_to_1");
		weightedOptions.bandwidth = 1.0;
	}
	if(weightedOptions.min_weight < 0.0)
	{
		weightedOptions.min_weight = 0.0;
	}
	append_hosc_window_tracef(
		"get_Weighted_Moment_Options.end enabled=%d mode=%d bandwidth=%.17g reference=%.17g min=%.17g normalize=%d zero_first=%d",
		useWeightedMoment ? 1 : 0,
		(int)weightedOptions.mode,
		weightedOptions.bandwidth,
		weightedOptions.reference_value,
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
	invertWindowOutput = parse_bool_text(invert_output, false);
	string enable_clip = parameter_value_compat(parameters, "enable_percentile_clipping.value");
	string clip_min = parameter_value_compat(parameters, "percentile_vmin.value");
	string clip_max = parameter_value_compat(parameters, "percentile_vmax.value");
	enablePercentileClipping = parse_bool_text(enable_clip, false);
	percentileClipMin = clip_min.empty() ? 2.0 : atof(clip_min.c_str());
	percentileClipMax = clip_max.empty() ? 98.0 : atof(clip_max.c_str());
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

	int iwl = atoi(wl.c_str());
	int iww = atoi(ww.c_str());
	int iwh = atoi(wh.c_str());

	int io = atoi(order.c_str());

	stepX = atoi(sx.c_str());
	stepY = atoi(sy.c_str());
	stepZ = atoi(sz.c_str());

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


	if(source_nx_ == 1 || source_ny_ == 1 || source_nz_ == 1)
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
	if(momentOrder < 0)
	{
		error_mesgs->report("moment_order","No less than 0");
		append_hosc_window_tracef("get_Windows_Moment.fail invalid_order order=%d", momentOrder);
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









	append_hosc_window_tracef(
		"get_Windows_Moment.end target_prop=%s parsed_window=%d,%d,%d parsed_step=%d,%d,%d output=%d,%d,%d invert=%d clip=%d clip_percentiles=%.17g,%.17g",
		auto_prop_name_.c_str(), momentWindowL, momentWindowW, momentWindowH, stepX, stepY, stepZ, myWindowShape.nx, myWindowShape.ny, myWindowShape.nz, invertWindowOutput ? 1 : 0, enablePercentileClipping ? 1 : 0, percentileClipMin, percentileClipMax);
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
	output_values.assign((size_t)gridSize, 0.0);
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
								weights.push_back(wmc::compute_weight(weighted_sample, weightedOptions));
							}
						}
					}
				}
				if(!complete_data_window || samples.empty())
				{
					output_values[(size_t)block] = 0.0;
					block++;
					continue;
				}
				output_keep_mask_[(size_t)block] = 1;
				kept_count++;

				double value = 0.0;
				if(useWeightedMoment)
				{
					try
					{
						if(G_CalculateTarget == 3)
						{
							value = wmc::weighted_univariate_moment(
								&samples[0], &weights[0], (int)samples.size(), momentOrder,
								weightedOptions.normalize_weights, weightedOptions.zero_first_order_moment);
						}
						else if(G_CalculateTarget == 4)
						{
							value = wmc::weighted_univariate_cumulant(
								&samples[0], &weights[0], (int)samples.size(), momentOrder,
								weightedOptions.normalize_weights, weightedOptions.zero_first_order_moment);
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
					if(al == Bell)
					{
						value = calculate_window_bell_cumulant_from_samples(samples, momentOrder);
					}
					else
					{
						if(!calculate_window_avl_cumulant_from_samples(samples, momentOrder, &value))
						{
							append_hosc_window_tracef("calculateMomentWindow.avl_failed block=%d sample_count=%d", block, (int)samples.size());
							return -1;
						}
					}
				}
				else
				{
					append_hosc_window_tracef("calculateMomentWindow.unsupported_target target=%d", G_CalculateTarget);
					return -1;
				}

				if(invertWindowOutput)
				{
					value = -value;
				}

				if(block == 0)
				{
					append_hosc_window_tracef(
						"calculateMomentWindow.first_block start=%d,%d,%d size=%d,%d,%d samples=%d weights=%d value=%.17g",
						start_x, start_y, start_z, length, width, height, (int)samples.size(), (int)weights.size(), value);
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
			if(i < (int)output_keep_mask_.size() && output_keep_mask_[(size_t)i] != 0)
				auto_property_->set_value(output_values[(size_t)i], i);
			else
				auto_property_->set_not_informed(i);
		}
		if(enablePercentileClipping)
		{
			const string clipped_name = auto_prop_name_ + "_clipped";
			GsTLGridProperty* clipped_property = geostat_utils::add_property_to_grid(auto_grid_, clipped_name);
			if(clipped_property != NULL)
			{
				for(int i = 0; i < value_count; i++)
				{
					if(i < (int)output_keep_mask_.size() && output_keep_mask_[(size_t)i] != 0)
						clipped_property->set_value(clipped_values[(size_t)i], i);
					else
						clipped_property->set_not_informed(i);
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

	_mkdir("D:/Development/DLL/WithSGeMS/MoveWindow/log");
	std::ostringstream tmp_stream;
	tmp_stream << "D:/Development/DLL/WithSGeMS/MoveWindow/log/"
		<< unique_name << "_" << (long)time(NULL) << ".sgems";
	const std::string tmp_filename = tmp_stream.str();
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
