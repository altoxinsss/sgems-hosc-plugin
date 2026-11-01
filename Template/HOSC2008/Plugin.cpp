
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
#include <errno.h>
#include <limits.h>
#include <limits>
#include <new>
#include <sstream>

#include "Plugin.h"
#include "PCA.h"
#include "DebugLogger.h"
#include "../../shared/ValidationFeedback.h"

namespace {
std::string normalize_algorithm_text(const std::string& text)
{
    std::string normalized = text;
    for(size_t i = 0; i < normalized.size(); i++)
    {
        normalized[i] = (char)tolower((unsigned char)normalized[i]);
    }
    return normalized;
}

// The SGeMS error handler inserts diagnostics into XML attributes verbatim.
// Escape parameter text before including it in a field-specific diagnostic.
std::string template_diagnostic_value(const std::string& text)
{
    std::string escaped;
    for(size_t i = 0; i < text.size(); ++i)
    {
        switch(text[i])
        {
        case '&': escaped += "&amp;"; break;
        case '<': escaped += "&lt;"; break;
        case '>': escaped += "&gt;"; break;
        case '"': escaped += "&quot;"; break;
        case '\'': escaped += "&apos;"; break;
        default: escaped += text[i]; break;
        }
    }
    return escaped;
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

bool parse_bool_selection(const std::string& text)
{
    std::string normalized = normalize_algorithm_text(text);
    return normalized == "true" || normalized == "yes" || normalized == "1" || normalized == "on" ||
        normalized.find("enable") != std::string::npos;
}

bool is_template_float_representable(double value)
{
	const double limit = (double)std::numeric_limits<float>::max();
	return weighted_moment_core::is_finite_number(value) && value >= -limit && value <= limit;
}

// Result rows retain one lag value for every reordered direction. SGeMS
// Cartesian properties are x-fastest, so rows are placed by lag indices, not
// by the enumeration order used to construct the lag table. The layout is
// captured before the calculator releases the Key-owned direction arrays.
int template_result_row_to_node_id(
	const TemplateOutputLayout& layout,
	const RGrid* output_grid,
	const double* row)
{
	if(!layout.valid || output_grid == NULL || row == NULL ||
	   (layout.displayed_direction_count != 2 && layout.displayed_direction_count != 3))
	{
		return -1;
	}
	const int dimensions[3] = { output_grid->nx(), output_grid->ny(), output_grid->nz() };
	for(int axis = 0; axis < 3; ++axis)
	{
		if(dimensions[axis] != layout.dimensions[axis] || dimensions[axis] <= 0)
		{
			return -1;
		}
	}
	int lag_index[3] = { 0, 0, 0 };
	for(int axis = 0; axis < layout.displayed_direction_count; ++axis)
	{
		const int column = layout.row_column[axis];
		const int step = layout.lag_step[axis];
		const double lag_value = row[column];
		if(column < 0 || step <= 0 || !weighted_moment_core::is_finite_number(lag_value))
		{
			return -1;
		}
		if(lag_value < 0.0 || lag_value > (double)INT_MAX - 0.5)
		{
			return -1;
		}
		const int rounded_lag = (int)floor(lag_value + 0.5);
		if(fabs(lag_value - (double)rounded_lag) > 1.0e-9 ||
		   rounded_lag % step != 0)
		{
			return -1;
		}
		lag_index[axis] = rounded_lag / step;
		if(lag_index[axis] < 0 || lag_index[axis] >= dimensions[axis])
		{
			return -1;
		}
	}
	if(layout.displayed_direction_count == 2 && dimensions[2] != 1)
	{
		return -1;
	}
	return lag_index[0] + dimensions[0] * (lag_index[1] + dimensions[1] * lag_index[2]);
}

bool parse_strict_template_int(const std::string& text, int* value)
{
	char* end = NULL;
	long parsed = 0;
	if(value == NULL || text.empty())
	{
		return false;
	}
	errno = 0;
	parsed = strtol(text.c_str(), &end, 10);
	if(end == text.c_str() || errno == ERANGE || parsed < INT_MIN || parsed > INT_MAX)
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

bool parse_exact_template_directions(
	const std::string& text,
	int expected_count,
	int coordinate_dimensions,
	std::vector<ExactTemplateDirection>* directions)
{
	if(directions == NULL || expected_count <= 0 || (coordinate_dimensions != 2 && coordinate_dimensions != 3))
	{
		return false;
	}
	directions->clear();
	std::stringstream input(text);
	std::string token;
	while(std::getline(input, token, ';'))
	{
		for(size_t character = 0; character < token.size(); ++character)
		{
			if(token[character] == ',' || token[character] == '(' || token[character] == ')' || token[character] == '[' || token[character] == ']')
			{
				token[character] = ' ';
			}
		}
		std::stringstream coordinates(token);
		std::string part_x;
		std::string part_y;
		std::string part_z;
		std::string extra;
		ExactTemplateDirection direction;
		direction.x = 0;
		direction.y = 0;
		direction.z = 0;
		if(!(coordinates >> part_x >> part_y) ||
		   (coordinate_dimensions == 3 && !(coordinates >> part_z)) ||
		   (coordinates >> extra) ||
		   !parse_strict_template_int(part_x, &direction.x) ||
		   !parse_strict_template_int(part_y, &direction.y) ||
		   (coordinate_dimensions == 3 && !parse_strict_template_int(part_z, &direction.z)) ||
		   (direction.x == 0 && direction.y == 0 && direction.z == 0))
		{
			directions->clear();
			return false;
		}
		directions->push_back(direction);
	}
	return (int)directions->size() == expected_count;
}

weighted_moment_core::WeightMode parse_weight_mode_selection(const std::string& text)
{
    std::string normalized = normalize_algorithm_text(text);
    if(normalized.find("coordinate") != std::string::npos || normalized.find("spatial") != std::string::npos)
    {
        return weighted_moment_core::WeightMode_CoordinateGaussian;
    }
    if(normalized.find("value") != std::string::npos)
    {
        return weighted_moment_core::WeightMode_ValueGaussian;
    }
    return weighted_moment_core::WeightMode_Uniform;
}

enum SlicePlaneMode
{
    SlicePlane_XY = 0,
    SlicePlane_XZ = 1,
    SlicePlane_YZ = 2
};

int parse_slice_plane_selection(const std::string& text)
{
    std::string normalized = normalize_algorithm_text(text);
    if(normalized == "xy")
    {
        return SlicePlane_XY;
    }
    if(normalized == "xz")
    {
        return SlicePlane_XZ;
    }
    if(normalized == "yz")
    {
        return SlicePlane_YZ;
    }
    return -1;
}

const char* slice_plane_label(int plane_mode)
{
    switch(plane_mode)
    {
    case SlicePlane_XY: return "xy";
    case SlicePlane_XZ: return "xz";
    case SlicePlane_YZ: return "yz";
    default: return "unknown";
    }
}

const char* slice_fixed_axis_label(int plane_mode)
{
    switch(plane_mode)
    {
    case SlicePlane_XY: return "z";
    case SlicePlane_XZ: return "y";
    case SlicePlane_YZ: return "x";
    default: return "axis";
    }
}

std::string parameter_value_or_empty(const Parameters_handler* parameters, const char* name)
{
    if(parameters == NULL || name == NULL)
    {
        hosc_debug_log("parameter.read.fail", "parameters=%p name=%p", parameters, name);
        return "";
    }

    hosc_debug_log("parameter.read.begin", "key=%s", name);
    std::string value;
    try
    {
        value = parameters->value(name);
    }
    catch(...)
    {
        hosc_debug_log("parameter.read.fail", "exception while reading key=%s", name);
        return "";
    }
    hosc_debug_log("parameter.read.end", "key=%s value=%s", name, value.c_str());
    return value;
}

double parameter_double_or(const Parameters_handler* parameters, const char* name, double fallback)
{
    std::string value = parameter_value_or_empty(parameters, name);
    if(value.empty())
    {
        hosc_debug_log("parameter.double.default", "key=%s fallback=%.17g", name == NULL ? "(null)" : name, fallback);
        return fallback;
    }
    return atof(value.c_str());
}

bool safe_register_project_object(GsTL_project* project, const std::string& object_name, unsigned long* exception_code)
{
    if(exception_code != NULL)
    {
        *exception_code = 0;
    }
    if(project == NULL || object_name.empty())
    {
        return false;
    }
    try
    {
        project->new_object(object_name);
        return true;
    }
    catch(const std::exception&)
    {
        if(exception_code != NULL)
        {
            *exception_code = 1;
        }
        return false;
    }
    catch(...)
    {
        if(exception_code != NULL)
        {
            *exception_code = 1;
        }
        return false;
    }
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
        out = "HOSC_OUT";
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

std::string join_int_values(const int* values, int count)
{
    if(values == NULL)
    {
        return "na";
    }
    std::ostringstream oss;
    const int safe_count = std::max(0, std::min(count, 10));
    for(int i = 0; i < safe_count; i++)
    {
        if(i > 0)
        {
            oss << "x";
        }
        oss << values[i];
    }
    return oss.str();
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
    case 1: return "CUM";
    case 2: return "MOM";
    case 3: return "WMOM";
    case 4: return "WCUM";
    default: return "OUT";
    }
}

const char* algorithm_name_token(Algorithm algorithm)
{
    return algorithm == Bell ? "BELL" : "AVL";
}

const char* weight_name_token(weighted_moment_core::WeightMode mode)
{
    switch(mode)
    {
    case weighted_moment_core::WeightMode_CoordinateGaussian: return "WCG";
    case weighted_moment_core::WeightMode_ValueGaussian: return "WVG";
    case weighted_moment_core::WeightMode_Uniform:
    default:
        return "WU";
    }
}
}











HOSC_TEMPLATE::HOSC_TEMPLATE()
{
	hosc_install_exception_filter();
	hosc_debug_log("HOSC_TEMPLATE.ctor", "construct plugin instance build=template_slice_stats_20260526_2230");
	source_geostat_grid_ = 0;
    source_grid_ = 0;
    source_node_count_ = 0;
    source_nx_ = 0;
    source_ny_ = 0;
    source_nz_ = 0;
    source_x = 0;
	source_y = 0;
	source_z = 0;
	source_val = 0;
    target_grid_ = 0;

	grid = NULL;
	Key = NULL;
	al = AVL;
	weighted_enabled_ = false;
	template_sampling_exact_ = false;
	integer_direction_vectors_.clear();
	current_project_ = NULL;
	auto_grid_ = NULL;
	auto_grid_name_.clear();
	slice_grid_ = NULL;
	slice_property_ = NULL;
	slice_grid_name_.clear();
	output_layout_.valid = false;
	output_layout_.displayed_direction_count = 0;
	for(int axis = 0; axis < 3; ++axis)
	{
		output_layout_.row_column[axis] = -1;
		output_layout_.lag_step[axis] = 0;
		output_layout_.dimensions[axis] = 0;
	}
#if 0
	momentWindowL = 0;
	momentWindowW = 0;
	momentWindowH = 0;
	momentOrder = 0;
#endif

	G_CalculateTarget = -1;
	slice_output_enabled_ = false;
	slice_plane_mode_ = SlicePlane_XY;
	slice_position_ = 1;
	slice_source_grid_name_.clear();
	slice_source_property_name_.clear();




}


bool HOSC_TEMPLATE::initialize( const Parameters_handler* parameters,
							   Error_messages_handler* error_mesgs )
{
	hosc_validation::FeedbackGuard feedback(error_mesgs, "hosc_template_validation_summary");
	hosc_install_exception_filter();
	hosc_debug_log("initialize.begin", "start");

	if( !get_Output_Slice_Options(parameters, error_mesgs))
	{
		hosc_debug_log("initialize.fail", "get_Output_Slice_Options failed");
		return false;
	}
	if(slice_output_enabled_)
	{
		if(!get_Output_Slice_Target(parameters, error_mesgs))
		{
			hosc_debug_log("initialize.fail", "get_Output_Slice_Target failed");
			return false;
		}
		if( !error_mesgs->empty() )
		{
			hosc_debug_log("initialize.fail", "slice mode error_mesgs is not empty");
			return false;
		}
		hosc_debug_log("initialize.end", "slice-only mode");
		return true;
	}

	if(!get_Calculate_Target(parameters, error_mesgs))
	{
		return false;
	}
	hosc_debug_log("initialize.target", "target=%d", G_CalculateTarget);


	if(!get_Calculate_Algorithm(parameters, error_mesgs) ||
	   !get_Weighted_Options(parameters, error_mesgs))
	{
		return false;
	}
	hosc_debug_log(
		"initialize.options",
		"algorithm=%d weighted=%d",
		(int)al,
		weighted_enabled_ ? 1 : 0);


    if ( !get_Source_Data(parameters, error_mesgs) )
	{
		hosc_debug_log("initialize.fail", "get_Source_Data failed");
        return false;
	}
	if ( !get_Keys_Data(parameters, error_mesgs))
	{
		hosc_debug_log("initialize.fail", "get_Keys_Data failed");
		return false;
	}
	if(!get_Template_Sampling(parameters, error_mesgs))
	{
		hosc_debug_log("initialize.fail", "get_Template_Sampling failed");
		return false;
	}



	if( !get_My_Shape(parameters, error_mesgs))
	{
		hosc_debug_log("initialize.fail", "get_My_Shape failed");
		return false;
	}
















	if( !error_mesgs->empty() )
	{
		hosc_debug_log("initialize.fail", "error_mesgs is not empty");
		return false;
	}



	return true;
}

int HOSC_TEMPLATE::execute( GsTL_project* proj )
{
	hosc_install_exception_filter();
	hosc_debug_log("execute.begin", "proj=%p source_geostat_grid=%p source_rgrid=%p grid=%p key=%p", proj, source_geostat_grid_, source_grid_, grid, Key);


	current_project_ = proj;
	if(slice_output_enabled_)
	{
		hosc_debug_log("execute.slice_only.begin", "grid=%s property=%s", slice_source_grid_name_.c_str(), slice_source_property_name_.c_str());
		if(!create_3d_slice_output_from_selected_grid())
		{
			GsTLcerr << "Failed to create 2D slice output from selected 3D grid." << gstlIO::end;
			hosc_debug_log("execute.slice_only.fail", "create_3d_slice_output_from_selected_grid failed");
			return 1;
		}
		hosc_debug_log("execute.slice_only.end", "success");
		return 0;
	}
	if(source_geostat_grid_ == NULL)
	{
		hosc_debug_log("execute.fail", "source_geostat_grid_ is NULL");
		return 1;
	}
	if(source_grid_ != NULL)
	{
		source_grid_->set_level(1);
	}



	if(grid == NULL)
	{

		hosc_debug_log("execute.fail", "grid is NULL");
		return 1;
	}
	if(Key == NULL)
	{
		hosc_debug_log("execute.fail", "Key is NULL");
		return 1;
	}
#if 0
	current_project_ = proj;
	auto_grid_ = CreateMyGrid(&myGridShape, "myGrid");
	slice_grid_ = NULL;
	slice_property_ = NULL;
	slice_grid_name_.clear();
	int gridSize = auto_grid_->nxyz();

	if(auto_grid_->nz() > 1 && auto_grid_->nx() > 1 && auto_grid_->ny() > 1)
		G_TargetDimension = 3;
	else
	{
		G_TargetDimension = 2;
		if(auto_grid_->nz() == 1)
			gridSize = auto_grid_->nx() * auto_grid_->ny();
		else if (auto_grid_->ny() == 1)
			gridSize = auto_grid_->nx() * auto_grid_->nz();
		else
			gridSize = auto_grid_->ny() * auto_grid_->nz();
	}
#endif

	double** targetGrid;


	int totalDim = Key->NDIR;
	for(int i = 0; i < Key->NDIR; i++)
	{
		if(Key->NLAG[i] <= 1)
			totalDim--;
	}
	int avlCrossection[10];
	for(int i = 0; i < crossectionNum && i < 10; i++)
	{
		avlCrossection[i] = crossection[i];
	}
	hosc_debug_log(
		"execute.layout",
		"target=%d algorithm=%d weighted=%d key_dim=%d ndir=%d order=%d output_shape=%d,%d,%d cross_count=%d cross0_zero=%d avl_cross0_zero=%d",
		G_CalculateTarget,
		(int)al,
		weighted_enabled_ ? 1 : 0,
		Key->DIM,
		Key->NDIR,
		Key->ORDER,
		myGridShape.nx,
		myGridShape.ny,
		myGridShape.nz,
		crossectionNum,
		crossectionNum > 0 ? crossection[0] : -1,
		crossectionNum > 0 ? avlCrossection[0] : -1);

	if(G_CalculateTarget == 1)
	{

		current_project_ = proj;
		hosc_debug_log("execute.create_grid.begin", "cumulant");
		auto_grid_ = CreateMyGrid(&myGridShape, "myGrid");
		slice_grid_ = NULL;
		slice_property_ = NULL;
		slice_grid_name_.clear();
		if(auto_grid_ == NULL || auto_property_ == NULL)
		{
			GsTLcerr << "Failed to create HOSC output grid or output property." << gstlIO::end;
			hosc_debug_log("execute.create_grid.fail", "auto_grid=%p auto_property=%p", auto_grid_, auto_property_);
			return 1;
		}
		int gridSize = auto_grid_->nxyz();
		hosc_debug_log("execute.create_grid.end", "gridSize=%d auto_nxyz=%d", gridSize, auto_grid_->nxyz());



		if(auto_grid_->nz() > 1 && auto_grid_->nx() > 1 && auto_grid_->ny() > 1)
			G_TargetDimension = 3;
		else
		{
			G_TargetDimension = 2;
			if(auto_grid_->nz() == 1)
				gridSize = auto_grid_->nx() * auto_grid_->ny();
			else if (auto_grid_->ny() == 1)
				gridSize = auto_grid_->nx() * auto_grid_->nz();
			else
				gridSize = auto_grid_->ny() * auto_grid_->nz();
		}

		if(template_sampling_exact_)
		{
			hosc_debug_log("execute.calc.begin", "ExactRecCalculation4DLL exact unweighted cumulant gridSize=%d", gridSize);
			targetGrid = ExactRecCalculation4DLL(
				Key,
				grid,
				&gridSize,
				avlCrossection,
				crossectionNum,
				integer_direction_vectors_,
				al == Bell);
		}
		else if(weighted_enabled_)
		{
			if(al == AVL)
			{
				hosc_debug_log("execute.calc.begin", "WeightedRecCalculation4DLL2 AVL weighted cumulant gridSize=%d", gridSize);
				targetGrid = WeightedRecCalculation4DLL2(
					Key,
					grid,
					&gridSize,
					avlCrossection,
					crossectionNum,
					weighted_options_,
					false);
			}
			else
			{
				hosc_debug_log("execute.calc.begin", "WeightedMoment4Dll Bell weighted cumulant gridSize=%d", gridSize);
				targetGrid = WeightedMoment4Dll(
					Key,
					grid,
					&gridSize,
					avlCrossection,
					crossectionNum,
					weighted_options_,
					false);
			}
		}
		else if(al == AVL)
		{


			hosc_debug_log("execute.calc.begin", "RecCalculation4DLL2 AVL cumulant gridSize=%d", gridSize);
			targetGrid = RecCalculation4DLL2(Key, grid, &gridSize, avlCrossection, crossectionNum);
		}
		else
		{


			hosc_debug_log("execute.calc.begin", "RecCalculation4DLL3 Bell cumulant gridSize=%d", gridSize);
			targetGrid = RecCalculation4DLL3(Key, grid, &gridSize, avlCrossection, crossectionNum);
		}
		hosc_debug_log("execute.calc.end", "targetGrid=%p returnedGridSize=%d grid_member=%p", targetGrid, gridSize, grid);

		if(targetGrid == NULL)
		{
			GsTLcerr << "HOSC template cumulant calculation returned no output. Check source property, NLAG/LAG ranges, and cross-section settings." << gstlIO::end;
			hosc_debug_log("execute.calc.fail", "targetGrid is NULL");
			return 1;
		}








#if 0
		if(totalDim > G_TargetDimension)
		{
#if 0
			double** trueTarget = (double**)calloc(gridSize, sizeof(double*));
			if(trueTarget == NULL)
			{
				return -1;
			}
			for(int i = 0; i < gridSize; i++)
			{
				trueTarget[i] = (double*)calloc(G_TargetDimension+1, sizeof(double));
				if(trueTarget[i] == NULL)
				{
					return -1;
				}
			}
#endif
			int output_rows = auto_grid_->nxyz();
			int output_cols = G_TargetDimension + 1;
			double **trueTarget = create_output_matrix(output_rows, output_cols);

			PCA_transform_advanced(targetGrid, trueTarget, gridSize, Key->ORDER, G_TargetDimension, auto_grid_->nxyz(), 0);
			set_MyGrid_Data(trueTarget, gridSize, G_TargetDimension);
		}
		else
#endif



			hosc_debug_log("execute.write.begin", "gridSize=%d valuePos=%d", gridSize, Key->NDIR);
			if(set_MyGrid_Data(targetGrid, gridSize, Key->NDIR) != 0)
				return 1;
			if(!write_Template_Output_Validity(targetGrid, gridSize, Key->NDIR))
				return 1;
			if(template_sampling_exact_ && !write_Template_Exact_Support(targetGrid, gridSize, Key->NDIR))
				return 1;
			hosc_debug_log("execute.write.end", "done");
		if(!create_3d_slice_output())
		{
			GsTLcerr << "Failed to create 2D slice output from 3D template result." << gstlIO::end;
		hosc_debug_log("execute.slice.fail", "slice output creation failed");
		return 1;
			}
			hosc_debug_log("execute.slice.end", "slice output ready name=%s", auto_grid_name_.c_str());

	}
	else if(G_CalculateTarget == 2)
	{
		current_project_ = proj;
		hosc_debug_log("execute.create_grid.begin", "moment");
		auto_grid_ = CreateMyGrid(&myGridShape, "myGrid");
		if(auto_grid_ == NULL || auto_property_ == NULL)
		{
			GsTLcerr << "Failed to create HOSC output grid or output property." << gstlIO::end;
			hosc_debug_log("execute.create_grid.fail", "moment auto_grid=%p auto_property=%p", auto_grid_, auto_property_);
			return 1;
		}
		int gridSize = auto_grid_->nxyz();
		hosc_debug_log("execute.create_grid.end", "moment gridSize=%d auto_nxyz=%d", gridSize, auto_grid_->nxyz());

		if(auto_grid_->nz() > 1 && auto_grid_->nx() > 1 && auto_grid_->ny() > 1)
			G_TargetDimension = 3;
		else
		{
			G_TargetDimension = 2;
			if(auto_grid_->nz() == 1)
				gridSize = auto_grid_->nx() * auto_grid_->ny();
			else if (auto_grid_->ny() == 1)
				gridSize = auto_grid_->nx() * auto_grid_->nz();
			else
				gridSize = auto_grid_->ny() * auto_grid_->nz();
		}
		if(weighted_enabled_)
		{
			hosc_debug_log("execute.calc.begin", "WeightedMoment4Dll moment gridSize=%d", gridSize);
			targetGrid = WeightedMoment4Dll(
				Key,
				grid,
				&gridSize,
				avlCrossection,
				crossectionNum,
				weighted_options_,
				true);
		}
		else
		{
			hosc_debug_log("execute.calc.begin", "Moment4Dll gridSize=%d", gridSize);
			targetGrid = Moment4Dll(Key, grid, &gridSize, avlCrossection, crossectionNum);
		}
		hosc_debug_log("execute.calc.end", "moment targetGrid=%p returnedGridSize=%d grid_member=%p", targetGrid, gridSize, grid);
		if(targetGrid == NULL)
		{
			GsTLcerr << "HOSC template moment calculation returned no output. Check source property, NLAG/LAG ranges, and cross-section settings." << gstlIO::end;
			hosc_debug_log("execute.calc.fail", "moment targetGrid is NULL");
			return 1;
		}
#if 0
		if(totalDim > G_TargetDimension)
		{
			int output_rows = auto_grid_->nxyz();
			int output_cols = G_TargetDimension + 1;
			double **trueTarget = create_output_matrix(output_rows, output_cols);

			PCA_transform_advanced(targetGrid, trueTarget, gridSize, Key->ORDER, G_TargetDimension, auto_grid_->nxyz(), 0);
			set_MyGrid_Data(trueTarget, auto_grid_->nxyz(), G_TargetDimension);
		}
		else
#endif
			hosc_debug_log("execute.write.begin", "moment gridSize=%d valuePos=%d", gridSize, Key->NDIR);
			if(set_MyGrid_Data(targetGrid, gridSize, Key->NDIR) != 0)
				return 1;
			if(!write_Template_Output_Validity(targetGrid, gridSize, Key->NDIR))
				return 1;
			hosc_debug_log("execute.write.end", "moment done");
			if(!create_3d_slice_output())
			{
				GsTLcerr << "Failed to create 2D slice output from 3D template result." << gstlIO::end;
				hosc_debug_log("execute.slice.fail", "moment slice output creation failed");
				return 1;
			}
			hosc_debug_log("execute.slice.end", "moment slice output ready name=%s", auto_grid_name_.c_str());
	}
#if 0
	else if(G_CalculateTarget == 3 || G_CalculateTarget == 4)
	{
		current_project_ = proj;

		auto_grid_ = CreateMyGrid(&myWindowShape, "myWM");

		if(calculateMomentWindow() == -1)
		{

			return -1;
		}
	}
#endif
	hosc_debug_log("execute.end", "success");
    return 0;
}

bool HOSC_TEMPLATE::get_Source_Data(const Parameters_handler* parameters,
							       Error_messages_handler* error_mesgs )
{

	string grid_name = parameters->value( "Source_Grid.value" );
	hosc_debug_log("get_Source_Data.begin", "grid_name=%s", grid_name.c_str());

    if ( grid_name.empty() )
    {
        error_mesgs->report( "Source_Grid", "No source grid selected" );
		hosc_debug_log("get_Source_Data.fail", "empty source grid");
        return false;
    }




    SmartPtr<Named_interface> source_ni =
        Root::instance()->interface(gridModels_manager + "/" + grid_name);
    source_geostat_grid_ = dynamic_cast<Geostat_grid*>(source_ni.raw_ptr());
    source_grid_ = dynamic_cast<RGrid*>(source_ni.raw_ptr());

    if( !source_geostat_grid_ )
    {
        std::ostringstream error_stream;
        error_stream <<  grid_name <<  " is not a valid geostatistical grid object";
        error_mesgs->report( "Source_Grid", error_stream.str() );
		hosc_debug_log("get_Source_Data.fail", "source geostat grid lookup failed");
        return false;
    }
	if(source_grid_ == NULL)
	{
		error_mesgs->report("Source_Grid", "HOSC_TEMPLATE supports Cartesian grids only");
		hosc_debug_log("get_Source_Data.fail", "source grid is not Cartesian");
		return false;
	}
    source_node_count_ = source_geostat_grid_->size();
    source_nx_ = source_grid_ != NULL ? source_grid_->nx() : 0;
    source_ny_ = source_grid_ != NULL ? source_grid_->ny() : 0;
    source_nz_ = source_grid_ != NULL ? source_grid_->nz() : 0;
	hosc_debug_log(
		"get_Source_Data.grid",
		"nodes=%d is_regular=%d nx=%d ny=%d nz=%d classname=%s",
		source_node_count_,
		source_grid_ != NULL ? 1 : 0,
		source_nx_,
		source_ny_,
		source_nz_,
		source_geostat_grid_->classname().c_str());


	// Cartesian-grid geometry is defined by node index and grid dimensions.
	// Ordinary x/y/z properties are data columns and must not affect templates.
	source_x = NULL;
	source_y = NULL;
	source_z = NULL;
	string prop_name = parameters->value( "Source_Property.value" );
	hosc_debug_log(
		"get_Source_Data.props",
		"source_prop=%s has_x=%d has_y=%d has_z=%d",
		prop_name.c_str(),
		source_x != NULL ? 1 : 0,
		source_y != NULL ? 1 : 0,
		source_z != NULL ? 1 : 0);
	if ( prop_name.empty() )
	{
		error_mesgs->report( "Source_Property", "No source property specified" );
		hosc_debug_log("get_Source_Data.fail", "empty source property");
		return false;
	}
	source_val = source_geostat_grid_->select_property(prop_name);
	if(source_val == NULL)
	{
		error_mesgs->report("Source_Property", "Selected source property does not exist");
		hosc_debug_log("get_Source_Data.fail", "source property does not exist: %s", prop_name.c_str());
		return false;
	}
	int num_xyz = source_node_count_;
	if(num_xyz <= 0)
	{
		error_mesgs->report("Source_Grid", "Selected source grid has no nodes");
		hosc_debug_log("get_Source_Data.fail", "source node count is zero");
		return false;
	}

	int informed_count = 0;
	for(int node = 0; node < num_xyz; node++)
	{
		if(!source_val->is_informed(node))
		{
			continue;
		}
		const double value = (double)source_val->get_value(node);
		if(!weighted_moment_core::is_finite_number(value))
		{
			continue;
		}
		informed_count++;
	}
	if(informed_count == 0)
	{
		error_mesgs->report("Source_Property", "Selected source property has no informed finite values");
		hosc_debug_log("get_Source_Data.fail", "no informed finite source values");
		return false;
	}

	grid = new (std::nothrow) double*[1+informed_count];

	if(grid == NULL)
	{
		error_mesgs->report("Source_Grid", "calloc return null");
		hosc_debug_log("get_Source_Data.fail", "grid row pointer allocation failed");
		return false;
	}
	grid[0] = new (std::nothrow) double[2];
	if(grid[0] == NULL)
	{
		error_mesgs->report("Source_Grid", "calloc return null");
		hosc_debug_log("get_Source_Data.fail", "grid header allocation failed");
		return false;
	}
	grid[0][0] = (double)informed_count;
	grid[0][1] = 4;
	int row = 1;
	for(int zero_index = 0; zero_index < num_xyz; zero_index++)
	{
		if(!source_val->is_informed(zero_index))
		{
			continue;
		}
		const double source_value = (double)source_val->get_value(zero_index);
		if(!weighted_moment_core::is_finite_number(source_value))
		{
			continue;
		}
		grid[row] = new (std::nothrow) double[4];
		if(grid[row] == NULL)
		{
			hosc_debug_log("get_Source_Data.fail", "grid row allocation failed at row=%d", row);
			return false;
		}
		const int nx = source_nx_;
		const int ny = source_ny_;
		const int safe_nx = (nx > 0) ? nx : 1;
		const int safe_ny = (ny > 0) ? ny : 1;
		const int xy = safe_nx * safe_ny;
		// This plugin boundary is always zero-based; it must not depend on the
		// first informed source node that happens to be encountered.
		grid[row][0] = (double)(zero_index % safe_nx);
		grid[row][1] = (double)((zero_index / safe_nx) % safe_ny);
		grid[row][2] = (double)(xy > 0 ? zero_index / xy : 0);
		grid[row][3] = source_value;
		row++;


	}
	if(num_xyz > 0)
	{
		hosc_debug_log(
			"get_Source_Data.sample",
			"valid=%d first=(%.6f,%.6f,%.6f,%.6f) last=(%.6f,%.6f,%.6f,%.6f)",
			informed_count,
			grid[1][0], grid[1][1], grid[1][2], grid[1][3],
			grid[informed_count][0], grid[informed_count][1], grid[informed_count][2], grid[informed_count][3]);
	}





	hosc_debug_log("get_Source_Data.end", "nodes=%d", num_xyz);
    return true;
}


bool HOSC_TEMPLATE::get_Target_Data(const Parameters_handler* parameters,
							         Error_messages_handler* error_mesgs )
{

    string grid_name = parameters->value( "Target_Grid.value" );

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











	target_prop_name_ = parameters->value( "props.value" );
	if(target_prop_name_.empty())
	{
		error_mesgs->report("props","No property specified");
		return false;
	}

#if 0
	target_dim_num = parameters->value("dims.value");
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

bool HOSC_TEMPLATE::get_Keys_Data(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	hosc_debug_log("get_Keys_Data.begin", "start");
	Key = (Keys*)calloc(1, sizeof(Keys));
	if(Key == NULL)
	{
		hosc_debug_log("get_Keys_Data.fail", "Key allocation failed");
		return false;
	}








	Key->REG_IRREG = 1;




	// Use the Cartesian grid's declared dimensions.  Inferring dimensions from
	// the informed node subset would shift or shrink geometry at missing edges.
	const int inferred_nx = source_nx_;
	const int inferred_ny = source_ny_;
	const int inferred_nz = source_nz_;
	Key->NX = inferred_nx;




	Key->NY = inferred_ny;




	Key->NZ = inferred_nz;
	hosc_debug_log("get_Keys_Data.inferred_dims", "nx=%d ny=%d nz=%d", Key->NX, Key->NY, Key->NZ);

	if(Key->NZ <= 1)
		Key->DIM = 2;
	else
		Key->DIM = 3;

#if 0
	string DIMX_TEXT = parameters->value("DIMX.value");
	stringstream ss6(DIMX_TEXT);
	ss6 >> Key->DIMX;

	string DIMY_TEXT = parameters->value("DIMY.value");
	stringstream ss7(DIMY_TEXT);
	ss7 >> Key->DIMY;

	string DIMZ_TEXT = parameters->value("DIMZ.value");
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




















	const char* direction_fields[] = {
		"NOD", "NLAG", "LAG", "LTOL", "AZM", "ATOL", "BANDWH", "DIP", "DTOL", "BANDWD"
	};
	vector<string>texts;
	texts.push_back(parameters->value("NOD.value"));
	texts.push_back(parameters->value("NLAG.value"));
	texts.push_back(parameters->value("LAG.value"));
	texts.push_back(parameters->value("LTOL.value"));
	texts.push_back(parameters->value("AZM.value"));
	texts.push_back(parameters->value("ATOL.value"));
	texts.push_back(parameters->value("BANDWH.value"));
	texts.push_back(parameters->value("DIP.value"));
	texts.push_back(parameters->value("DTOL.value"));
	texts.push_back(parameters->value("BANDWD.value"));
	hosc_debug_log(
		"get_Keys_Data.raw",
		"NOD=%s NLAG=%s LAG=%s LTOL=%s AZM=%s",
		texts[0].c_str(),
		texts[1].c_str(),
		texts[2].c_str(),
		texts[3].c_str(),
		texts[4].c_str());

	int dirNum = 0;
	vector<string> tempResult;
	for(int i = 0; i < texts.size(); i++)
	{
		if(texts[i].empty())
		{
			error_mesgs->report(direction_fields[i], "Template shape: this direction parameter list is empty. Enter one value per direction, separated by commas.");
			return false;
		}
		stringstream ss(texts[i]);
		string strItem;
		int tempNum = 0;
		while(getline(ss, strItem, ','))
		{
			if(strItem.empty())
			{
				error_mesgs->report(direction_fields[i], "Template shape: this list contains an empty item. Enter one value per direction, separated by commas, without empty items.");
				return false;
			}
			tempResult.push_back(strItem);
			tempNum++;
			if(i == 0)
				dirNum++;
		}
		if(tempNum != dirNum)
		{
			std::ostringstream diagnostic;
			diagnostic << "Template shape: " << direction_fields[i] << " contains " << tempNum
				<< " values, but Number of occurrences (NOD) defines " << dirNum
				<< " directions. Enter exactly one value per direction.";
			error_mesgs->report(direction_fields[i], diagnostic.str());
			hosc_debug_log("get_Keys_Data.fail", "parameter count mismatch group=%d tempNum=%d dirNum=%d", i, tempNum, dirNum);
			return false;
		}
	}
	if(dirNum <= 0 || dirNum > 10)
	{
		error_mesgs->report("NOD", "The number of directions must be between 1 and 10");
		return false;
	}
	Key->NDIR = dirNum;
	Key->ORDER = dirNum+1;
	hosc_debug_log("get_Keys_Data.order", "ndir=%d order=%d", Key->NDIR, Key->ORDER);




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

	for(int group = 0; group < (int)tempInt.size(); group++)
	{
		for(int j = 0; j < Key->NDIR; j++)
		{
			const int value_index = group * Key->NDIR + j;
			int parsed = 0;
			if(value_index >= (int)tempResult.size() || !parse_strict_template_int(tempResult[value_index], &parsed))
			{
				std::ostringstream diagnostic;
				diagnostic << "Template shape: " << direction_fields[group] << ", direction " << (j + 1)
					<< " requires a complete integer without trailing characters.";
				error_mesgs->report(direction_fields[group], diagnostic.str());
				return false;
			}
			tempInt[group][j] = parsed;
		}
	}
	bool valid_direction_ranges = true;
	for(int group = 0; group < (int)tempInt.size(); group++)
	{
		// Azimuth and dip may be signed; only counts, steps, and tolerances
		// have the positivity constraints already enforced by this parser.
		if(group == 4 || group == 7)
		{
			continue;
		}
		for(int direction = 0; direction < Key->NDIR; direction++)
		{
			const int value = tempInt[group][direction];
			if((group < 3 && value <= 0) || (group >= 3 && value < 0))
			{
				std::ostringstream diagnostic;
				diagnostic << "Template shape: " << direction_fields[group] << ", direction " << (direction + 1)
					<< " has value " << value << ". Enter "
					<< (group < 3 ? "a positive integer." : "a non-negative integer.");
				error_mesgs->report(direction_fields[group], diagnostic.str());
				valid_direction_ranges = false;
				break;
			}
		}
	}
	if(!valid_direction_ranges)
	{
		return false;
	}
	int requested_joint_order = 1;
	for(int direction = 0; direction < Key->NDIR; direction++)
	{
		requested_joint_order += Key->NOD[direction];
	}
	if(requested_joint_order > 8)
	{
		error_mesgs->report("NOD", "The supported total cumulant order is at most 8");
		return false;
	}
	vector<string>().swap(texts);
	vector<string>().swap(tempResult);
	if(Key->NDIR > 0)
	{
		hosc_debug_log(
			"get_Keys_Data.end",
			"first_dir nod=%d nlag=%d lag=%d azm=%d dims=%d,%d,%d dim=%d",
			Key->NOD[0],
			Key->NLAG[0],
			Key->LAG[0],
			Key->AZM[0],
			Key->NX,
			Key->NY,
			Key->NZ,
			Key->DIM);
	}









	return true;
}


int HOSC_TEMPLATE::set_Target_Data(double** targetGrid, int gridSize, int valuePos)
{

	for(int i = 0; i < gridSize; i++)
	{
		target_property_->set_value(targetGrid[i][valuePos], i);

	}


	return 0;
}


bool HOSC_TEMPLATE::get_Calculate_Target(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	string Calculate_Target = parameter_value_or_empty(parameters, "Target_Button.value");

	if(Calculate_Target == "Moment")
		G_CalculateTarget = 2;
	else if(Calculate_Target == "Cumulant")
		G_CalculateTarget = 1;
	else
	{
		error_mesgs->report("Target_Button", "Expected Moment or Cumulant");
		return false;
	}
	hosc_debug_log("get_Calculate_Target", "raw=%s parsed=%d", Calculate_Target.c_str(), G_CalculateTarget);

	return true;
}

bool HOSC_TEMPLATE::get_Calculate_Algorithm(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	string Calculate_Algorithm = parameter_value_or_empty(parameters, "Algorithm_Mode.value");
	if(Calculate_Algorithm.empty())
	{
		Calculate_Algorithm = parameter_value_or_empty(parameters, "Algorithm_Button.value");
	}
	const std::string normalized = normalize_algorithm_text(Calculate_Algorithm);
	if(normalized == "bell" || normalized == "bell algorithm" || normalized == "bell based algorithm")
	{
		al = Bell;
	}
	else if(normalized == "avl" || normalized == "avl algorithm" || normalized == "avl based algorithm")
	{
		al = AVL;
	}
	else
	{
		error_mesgs->report("Algorithm_Mode", "General / Calculate Algorithm: select Bell or AVL.");
		return false;
	}
	hosc_debug_log("get_Calculate_Algorithm", "raw=%s parsed=%d", Calculate_Algorithm.c_str(), (int)al);
	return true;
}

bool HOSC_TEMPLATE::get_Template_Sampling(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	template_sampling_exact_ = false;
	integer_direction_vectors_.clear();
	const std::string sampling_text = parameter_value_or_empty(parameters, "template_sampling.value");
	const std::string sampling = normalize_algorithm_text(sampling_text);
	const bool legacy_sampling = sampling.empty() || sampling == "legacy" || sampling == "legacy_search";
	if(legacy_sampling && source_nz_ <= 1)
	{
		return true;
	}
	if(!legacy_sampling && sampling != "exact" && sampling != "exhaustive_integer_offsets")
	{
		std::ostringstream diagnostic;
		diagnostic << "Exact Template / Sampling mode: current value '" << template_diagnostic_value(sampling_text)
			<< "' is not supported. Select legacy_search for a 2D grid, or exhaustive_integer_offsets for an unweighted Cumulant.";
		error_mesgs->report("template_sampling", diagnostic.str());
		return false;
	}

	// A 3D grid requires exact sampling. Report all incompatible inputs at
	// once so correcting the mode does not merely reveal the next error.
	bool valid = true;
	if(legacy_sampling)
	{
		std::ostringstream diagnostic;
		diagnostic << "Exact Template / Sampling mode: current value '"
			<< template_diagnostic_value(sampling_text.empty() ? "legacy_search (default)" : sampling_text)
			<< "' cannot sample this 3D Cartesian grid (Nz=" << source_nz_
			<< "). Select exhaustive_integer_offsets, then satisfy the other fields listed below.";
		error_mesgs->report("template_sampling", diagnostic.str());
		valid = false;
	}
	if(G_CalculateTarget != 1)
	{
		std::ostringstream diagnostic;
		diagnostic << "General / Calculate Target: current value '"
			<< template_diagnostic_value(parameter_value_or_empty(parameters, "Target_Button.value"))
			<< "' is incompatible with exhaustive_integer_offsets. Select Cumulant; exact template sampling does not support Moment output.";
		error_mesgs->report("Target_Button", diagnostic.str());
		valid = false;
	}
	if(weighted_enabled_)
	{
		error_mesgs->report("Weighted_Enable", "Template calculations are unweighted. Remove Weighted_Enable from the parameter file or set it to Disable.");
		valid = false;
	}
	const int coordinate_dimensions = source_nz_ <= 1 ? 2 : 3;
	const int direction_count = Key != NULL ? Key->NDIR : 0;
	const std::string direction_text = parameter_value_or_empty(parameters, "integer_direction_vectors.value");
	std::vector<ExactTemplateDirection> parsed_directions;
	if(!parse_exact_template_directions(direction_text, direction_count, coordinate_dimensions, &parsed_directions))
	{
		std::ostringstream diagnostic;
		diagnostic << "Exact Template / Integer direction vectors: current value '" << template_diagnostic_value(direction_text)
			<< "' is invalid. Enter exactly " << direction_count << " nonzero integer vectors, one per Template shape direction, with "
			<< coordinate_dimensions << " components each (" << (coordinate_dimensions == 3 ? "dx,dy,dz" : "dx,dy")
			<< "), separated by semicolons. Use directions appropriate to the intended template.";
		error_mesgs->report("integer_direction_vectors", diagnostic.str());
		valid = false;
	}
	if(!valid)
	{
		return false;
	}
	integer_direction_vectors_.swap(parsed_directions);
	template_sampling_exact_ = true;
	hosc_debug_log("get_Template_Sampling", "exhaustive integer offsets directions=%d", (int)integer_direction_vectors_.size());
	return true;
}

bool HOSC_TEMPLATE::get_Weighted_Options(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	hosc_debug_log("get_Weighted_Options.begin", "start");
	// All public template calculations, including legacy 2D sampling, use
	// unweighted moments. Validate old PAR files before any sampling branch.
	weighted_enabled_ = false;
	weighted_options_ = weighted_moment_core::WeightedMomentOptions();
	weighted_options_.normalize_weights = true;
	weighted_options_.zero_first_order_moment = false;
	const string enabled_text = normalize_algorithm_text(parameter_value_or_empty(parameters, "Weighted_Enable.value"));
	if(!enabled_text.empty() && enabled_text != "disable" && enabled_text != "disabled" &&
	   enabled_text != "false" && enabled_text != "0" && enabled_text != "no" && enabled_text != "off")
	{
		error_mesgs->report("Weighted_Enable", "Template calculations are unweighted. Remove Weighted_Enable from the parameter file or set it to Disable. Use HOSC_WINDOW for normalized spatial or value weighting.");
		hosc_debug_log("get_Weighted_Options.fail", "template weighting request rejected");
		return false;
	}
	hosc_debug_log("get_Weighted_Options.end", "fixed unweighted policy");
	return true;
}

#if 0
bool HOSC::get_Dim_Data(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	string Dim = parameters->value("dims.value");
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


#if 0
int HOSC_TEMPLATE::calculateMomentWindow()
{
	if(auto_grid_ == NULL)
	{

		return -1;
	}
	if(momentWindowL == 0 && momentWindowW == 0 && momentWindowH == 0)
	{

		return -1;
	}

	int gridSize = auto_grid_->nxyz();

	if(auto_grid_->nz() > 1 && auto_grid_->nx() > 1 && auto_grid_->ny() > 1)
		G_TargetDimension = 3;
	else
	{
		G_TargetDimension = 2;

	}



	if(Key->DIM == 2)
	{



		if(stepX == -1 && stepY == -1)
		{
			int i = 0;
			int start = 1;
			while(i < gridSize)
			{
				int lenght = momentWindowL;
				int width = momentWindowW;
				int height = momentWindowH;

				if(start > source_node_count_)



				if(grid[start][0] == 1)
				{
					lenght /= 2;
				}
				else if(grid[start][0] + lenght - 1 > Key->NX)
				{
					lenght = Key->NX - grid[start][0] + 1;
				}

				if(grid[start][1] == 1)
				{
					width /= 2;
				}
				else if(grid[start][1] + width - 1 > Key->NY)
				{
					width = Key->NY - grid[start][1] + 1;
				}

				int totalNode = lenght * width;


				int tempStart = start;
				double tempSum = 0;
				for(int w = 0; w < width; w++)
				{






					for(int l = 0; l < lenght; l++)
					{

						if(tempStart > source_node_count_)
						{



						}

						tempSum += pow(grid[tempStart][3], (double)momentOrder);
						tempStart++;
					}
					tempStart = tempStart + Key->NX - lenght;
				}
				tempSum /= totalNode;
				auto_property_->set_value(tempSum, i);
				i++;
				if(i % auto_grid_->nx() == 0)
				{

					start = start + lenght + (width-1) * Key->NX;



				}
				else
				{
					start += lenght;
				}
			}
		}
		else
		{
			int block = 1;

			int ld = 1;
			int startX = 0;
			int startY = 0;
			while(block <= gridSize)
			{
				int lenght = momentWindowL;
				int width = momentWindowW;
				int height = momentWindowH;
				startX = (grid[middle][0] - (lenght / 2) < 1) ? (1 : grid[middle][0] - (lenght / 2));
				startY = (grid[middle][1] - (width / 2) < 1) ? (1 : grid[middle][1] - (width / 2));
				lenght = grid[middle][0] + lenght / 2 - startX + 1;
			}
		}
	}

	return 0;
}

#endif


#if 0
bool HOSC_TEMPLATE::get_Windows_Moment(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{


	string wl = parameters->value("windowl.value");
	string ww = parameters->value("windoww.value");
	string wh = parameters->value("windowh.value");

	string order = parameters->value("moment_order.value");

	string sx = parameters->value("stepx.value");
	string sy = parameters->value("stepy.value");
	string sz = parameters->value("stepz.value");




	auto_prop_name_ = parameters->value( "props.value" );
	if(auto_prop_name_.empty())
	{
		error_mesgs->report("props","No property specified");
		return false;
	}

	int iwl = atoi(wl.c_str());
	int iww = atoi(ww.c_str());
	int iwh = atoi(wh.c_str());

	int io = atoi(order.c_str());

	stepX = atoi(sx.c_str());
	stepY = atoi(sy.c_str());
	stepZ = atoi(sz.c_str());


	if(Key->DIM == 2)
	{
		if((iwh != 1 && iwh != 0))
		{
			error_mesgs->report("windowz","should be 1 if souce grid is 2D");
			return false;
		}
		if(stepZ != 0 && stepZ != -1)
		{
			error_mesgs->report("stepz","should be 0 if souce grid is 2D");
			return false;
		}
	}



	if(iwh == 0)
	{
		iwh = 1;
		stepZ = -1;
	}

	momentWindowL = iwl;
	momentWindowW = iww;
	momentWindowH = iwh;
	momentOrder = io;
	if(momentOrder < 0)
	{
		error_mesgs->report("moment_order","No less than 0");
		return false;
	}

	int iinx = 0;
	int iiny = 0;
	int iinz = 0;

	if(stepX == -1)
	{
		iinx = (Key->NX - (momentWindowL / 2)) / momentWindowL + 2;
		if(Key->NX <= (momentWindowL / 2))
			iinx = 1;
		else if(Key->NX <= momentWindowL)
			iinx = 2;
		else if((Key->NX - (momentWindowL / 2)) % momentWindowL == 0)
			iinx--;
	}
	else if(stepX <= 0)
	{
		error_mesgs->report("stepx", "No less than 0");
		return false;
	}
	else
	{

		iinx = (Key->NX - 1) / stepX;
	}


	if(stepY == -1)
	{
		iiny = (Key->NY - (momentWindowW / 2)) / momentWindowW + 2;
		if(Key->NY <= (momentWindowW / 2))
			iiny = 1;
		else if(Key->NY <= momentWindowW)
			iiny = 2;
		else if((Key->NY - (momentWindowW / 2)) % momentWindowW == 0)
			iiny--;
	}
	else if(stepY <= 0)
	{
		error_mesgs->report("stepy", "No less than 0");
		return false;
	}
	else
	{

		iiny = (Key->NY - 1) / stepY;
	}


	if(stepZ == -1)
	{
		iinz = (Key->NZ - (momentWindowH / 2)) / momentWindowH + 2;
		if(Key->NZ <= (momentWindowH / 2))
			iinz = 1;
		else if(Key->NZ <= momentWindowH)
			iinz = 2;
		else if((Key->NZ - (momentWindowH / 2)) % momentWindowH == 0)
			iinz--;
	}
	else if(stepZ <= 0)
	{
		error_mesgs->report("stepz", "No less than 0");
		return false;
	}
	else
	{

		iinz = (Key->NZ - 1) / stepZ;
	}



	if(iinx <= 0 || iiny <= 0 || iinz <= 0)
	{
		error_mesgs->report("windowl","No less than 0");
		return false;
	}

	GridShape tempShape = {iinx, iiny, iinz, 0, 0, 0, momentWindowL, momentWindowW, momentWindowH};


	myWindowShape = tempShape;



	return true;
}
#endif

#if 0
bool HOSC_TEMPLATE::get_Windows_Moment(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	GsTLlog << "Start get_Windows_Moment" << gstlIO::end;

	string wl = parameters->value("windowl.value");
	string ww = parameters->value("windoww.value");
	string wh = parameters->value("windowh.value");

	string order = parameters->value("moment_order.value");

	string sx = parameters->value("stepx.value");
	string sy = parameters->value("stepy.value");
	string sz = parameters->value("stepz.value");




	auto_prop_name_ = parameters->value( "props.value" );
	if(auto_prop_name_.empty())
	{
		error_mesgs->report("props","No property specified");
		return false;
	}

	int iwl = atoi(wl.c_str());
	int iww = atoi(ww.c_str());
	int iwh = atoi(wh.c_str());

	int io = atoi(order.c_str());

	stepX = atoi(sx.c_str());
	stepY = atoi(sy.c_str());
	stepZ = atoi(sz.c_str());

	if(stepX <= 0 || stepY <= 0)
	{
		error_mesgs->report("stepx","No less than 0");
		return false;
	}
	if(Key->DIM == 3 && stepZ <= 0)
	{
		error_mesgs->report("stepz", "No less than 0");
	}


	if(Key->DIM == 2)
	{
		if((iwh != 1 && iwh != 0))
		{
			error_mesgs->report("windowz","should be 1 if souce grid is 2D");
			return false;
		}
		if(stepZ != 0)
		{
			error_mesgs->report("stepz","should be 0 if souce grid is 2D");
			return false;
		}
		iwh = 1;
	}


	momentWindowL = iwl;
	momentWindowW = iww;
	momentWindowH = iwh;
	momentOrder = io;
	if(momentOrder < 0)
	{
		error_mesgs->report("moment_order","No less than 0");
		return false;
	}

	int iinx = 0;
	int iiny = 0;
	int iinz = 0;


	iinx = 1 + (Key->NX - 1) / stepX;
	iiny = 1 + (Key->NY - 1) / stepY;
	if(Key->DIM == 2)
		iinz = 1;
	else
		iinz = 1 + (Key->NZ - 1) / stepZ;

	if(iinx <= 0 || iiny <= 0)
	{
		error_mesgs->report("windowl","No less than 0");
		return false;
	}
	if(Key->DIM == 3 && iinz <= 0)
	{
		error_mesgs->report("windowH", "No less than 0");
	}

	GridShape tempShape = {iinx, iiny, iinz, 0, 0, 0, momentWindowL, momentWindowW, momentWindowH};


	myWindowShape = tempShape;



	return true;
}


int HOSC_TEMPLATE::calculateMomentWindow()
{

	if(auto_grid_ == NULL)
	{

		return -1;
	}
	if(momentWindowL == 0 && momentWindowW == 0 && momentWindowH == 0)
	{

		return -1;
	}

	int gridSize = auto_grid_->nxyz();

	if(auto_grid_->nz() > 1 && auto_grid_->nx() > 1 && auto_grid_->ny() > 1)
		G_TargetDimension = 3;
	else
	{
		G_TargetDimension = 2;

	}




	int block = 0;
	int startX = 1;
	int startY = 1;
	int startZ = 1;

	while(block < gridSize)
	{
		int lenght = momentWindowL;
		int width = momentWindowW;
		int height = momentWindowH;

		int start = (startZ - 1) * Key->NX *Key->NY + (startY - 1) * Key->NX + startX;

		if(startX + lenght - 1 > Key->NX)
		{
			lenght = Key->NX - startX + 1;
		}
		if(startY + width - 1 > Key->NY)
		{
			width = Key->NY - startY + 1;
		}
		if(Key->DIM == 3 && startZ + height - 1 > Key->NZ)
		{
			height = Key->NZ - startZ + 1;
		}

		int totalNode = lenght * width * height;



		if(start > grid[0][0])
		{

			return -1;
		}

		if(G_CalculateTarget == 3)
		{
			int currNode = start;
			int currPos = 1;
			double tempSum = 0;

			for(int w = 1; w <= totalNode; w++)
			{
				tempSum += pow(grid[currNode][3], (double)momentOrder);

				currNode++;
				currPos++;
				if(currPos == lenght+1)
				{
					currNode = currNode + Key->NX - lenght;
					currPos = 1;
				}
			}
			tempSum /= totalNode;
			auto_property_->set_value(tempSum, block);
		}
		else if(G_CalculateTarget == 4)
		{

		}

		block++;
		startX += stepX;;
		if(startX > Key->NX)
		{
			startX = 1;
			startY += stepY;
			if(startY > Key->NY)
			{
				if(Key->DIM == 3)
				{
					startY = 1;
					startZ += stepZ;
					if(startZ > Key->NZ)
					{

						return 0;
					}
				}
				else if(Key->DIM == 2)
				{

					return 0;
				}
			}
		}
	}
	return 0;
}
#endif


bool HOSC_TEMPLATE::get_Sliding_Window(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	string wx = parameters->value("windowx.value");
	string wy = parameters->value("windowy.value");
	string wz = parameters->value("windowz.value");
	string wl = parameters->value("windowl.value");
	string ww = parameters->value("windoww.value");
	string wh = parameters->value("windowh.value");

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
		if(currNode >= 1+source_node_count_)
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





	Free2DVec(&grid, 1+source_node_count_);
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


#if 0
RGrid* HOSC_TEMPLATE::CreateMyGrid(GridShape* tempGridShape, const char* gridName)
{

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


    SmartPtr<Named_interface> ni =
        Root::instance()->new_interface("cgrid", full_path);

    if (ni.raw_ptr() == 0) {

        return NULL;
    }

    Cartesian_grid* grid = dynamic_cast<Cartesian_grid*>(ni.raw_ptr());

    if (!grid) {

        return NULL;
    }

    grid->set_dimensions(tempGridShape->nx, tempGridShape->ny, tempGridShape->nz,
                        tempGridShape->dx, tempGridShape->dy, tempGridShape->dz);

    grid->origin(GsTLPoint(tempGridShape->cx, tempGridShape->cy, tempGridShape->cz));


    if (current_project_) {
        current_project_->new_object(unique_name);

    }



    return grid;
}
#endif

void HOSC_TEMPLATE::TestMyGrid()
{
	GridShape myGridShape = {10, 10, 10, 0, 0, 0, 1, 1, 1};
	RGrid* myGrid = CreateMyGrid(&myGridShape, "MyGrid");
}

#if 0
RGrid* HOSC_TEMPLATE::CreateMyGrid(GridShape* tempGridShape, const char* gridName)
{
	if(Key == NULL)
	{

		return NULL;
	}



#if 0
	string scx = parameters->value( "xmn_.value" );
	string scy = parameters->value( "ymn_.value" );
	string scz = parameters->value( "zmn_.value" );
	string sdx = parameters->value( "x_size_" );
	string sdy = parameters->value( "y_size_" );
	string sdz = parameters->value( "z_size_" );

	double dcx = String_Op::to_number<double>(scx);
	double dcy = String_Op::to_number<double>(scy);
	double dcz = String_Op::to_number<double>(scz);
	double ddx = String_Op::to_number<double>(sdx);
	double ddy = String_Op::to_number<double>(sdy);
	double ddz = String_Op::to_number<double>(sdz);
#endif





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


    SmartPtr<Named_interface> ni =
        Root::instance()->new_interface("cgrid", full_path);

    if (ni.raw_ptr() == 0) {

        return NULL;
    }

    Cartesian_grid* myGrid = dynamic_cast<Cartesian_grid*>(ni.raw_ptr());

    if (!myGrid) {

        return NULL;
    }



    myGrid->set_dimensions(tempGridShape->nx, tempGridShape->ny, tempGridShape->nz,
                        tempGridShape->dx, tempGridShape->dy, tempGridShape->dz);



    myGrid->origin(GsTLPoint(tempGridShape->cx, tempGridShape->cy, tempGridShape->cz));









#if 0
    if (current_project_) {

        current_project_->new_object(unique_name);

    }
	else
	{

		return NULL;
	}
#endif





    if (current_project_) {
        try {







            current_project_->new_object(unique_name);


        } catch (const std::exception& e) {

            return NULL;
        } catch (...) {

            return NULL;
        }
    } else {

        return NULL;
    }






	auto_property_ = geostat_utils::add_property_to_grid(myGrid, auto_prop_name_);



	myGrid->set_level(1);

    return myGrid;
}
#endif


RGrid* HOSC_TEMPLATE::CreateMyGrid(GridShape* tempGridShape, const char* gridName)
{
	hosc_debug_log(
		"CreateMyGrid.begin",
		"gridName=%s shape=%d,%d,%d origin=%.3f,%.3f,%.3f cell=%.3f,%.3f,%.3f project=%p",
		gridName == NULL ? "(null)" : gridName,
		tempGridShape != NULL ? tempGridShape->nx : -1,
		tempGridShape != NULL ? tempGridShape->ny : -1,
		tempGridShape != NULL ? tempGridShape->nz : -1,
		tempGridShape != NULL ? tempGridShape->cx : 0.0,
		tempGridShape != NULL ? tempGridShape->cy : 0.0,
		tempGridShape != NULL ? tempGridShape->cz : 0.0,
		tempGridShape != NULL ? tempGridShape->dx : 0.0,
		tempGridShape != NULL ? tempGridShape->dy : 0.0,
		tempGridShape != NULL ? tempGridShape->dz : 0.0,
		current_project_);
    if(Key == NULL) {

		hosc_debug_log("CreateMyGrid.fail", "Key is NULL");
        return NULL;
    }



    uintptr_t proj_addr = reinterpret_cast<uintptr_t>(current_project_);


    if (proj_addr < 0x1000) {


		hosc_debug_log("CreateMyGrid.fail", "invalid current_project address=%p", current_project_);
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

		hosc_debug_log("CreateMyGrid.fail", "failed to generate unique name");
        return NULL;
    }
	hosc_debug_log("CreateMyGrid.name", "unique=%s path=%s", unique_name.c_str(), full_path.c_str());


    SmartPtr<Named_interface> ni =
        Root::instance()->new_interface("cgrid", full_path);

    if (ni.raw_ptr() == 0) {

		hosc_debug_log("CreateMyGrid.fail", "Root::new_interface returned NULL");
        return NULL;
    }

    Cartesian_grid* myGrid = dynamic_cast<Cartesian_grid*>(ni.raw_ptr());

    if (!myGrid) {

		hosc_debug_log("CreateMyGrid.fail", "dynamic_cast<Cartesian_grid*> failed");
        return NULL;
    }

    myGrid->set_dimensions(tempGridShape->nx, tempGridShape->ny, tempGridShape->nz,
                          tempGridShape->dx, tempGridShape->dy, tempGridShape->dz);

    myGrid->origin(GsTLPoint(tempGridShape->cx, tempGridShape->cy, tempGridShape->cz));







    try {

        GsTLlog << "Registering object: " << unique_name << gstlIO::end;
		hosc_debug_log("CreateMyGrid.register.begin", "unique=%s", unique_name.c_str());
        current_project_->new_object(unique_name);
        GsTLlog << "Grid registered successfully" << gstlIO::end;
		hosc_debug_log("CreateMyGrid.register.end", "registered");
    } catch (const std::exception& e) {

		hosc_debug_log("CreateMyGrid.fail", "std exception in new_object: %s", e.what());
        return NULL;
    } catch (...) {

		hosc_debug_log("CreateMyGrid.fail", "unknown exception in new_object");
        return NULL;
    }


    SmartPtr<Named_interface> managed_ni = Root::instance()->interface(full_path);
    if (!managed_ni.raw_ptr()) {

		hosc_debug_log("CreateMyGrid.fail", "failed to retrieve registered grid path=%s", full_path.c_str());
        return NULL;
    }



	RGrid* managedGrid = dynamic_cast<RGrid*>(managed_ni.raw_ptr());
	if(managedGrid == NULL)
	{
		hosc_debug_log("CreateMyGrid.fail", "dynamic_cast<RGrid*> failed");
		return NULL;
	}
	auto_property_ = geostat_utils::add_property_to_grid(managedGrid, auto_prop_name_);
	if(auto_property_ == NULL)
	{
		hosc_debug_log("CreateMyGrid.fail", "add_property_to_grid failed prop=%s", auto_prop_name_.c_str());
		return NULL;
	}
	managedGrid->set_level(1);
	auto_grid_name_ = unique_name;
	hosc_debug_log("CreateMyGrid.end", "managedGrid=%p property=%p nxyz=%d", managedGrid, auto_property_, managedGrid->nxyz());

    return managedGrid;
}

std::string HOSC_TEMPLATE::build_output_grid_base_name(const GridShape* tempGridShape, const char* role) const
{
    const std::string role_text = role == NULL ? "" : normalize_algorithm_text(role);
    const bool window_role = role_text.find("window") != std::string::npos;

    const char target_token =
        (G_CalculateTarget == 1 || G_CalculateTarget == 4) ? 'C' :
        (G_CalculateTarget == 2 || G_CalculateTarget == 3) ? 'M' : 'O';
    const char algorithm_token = al == Bell ? 'B' : 'A';
    char weight_token = 'U';
    if(weighted_enabled_)
    {
        if(weighted_options_.mode == weighted_moment_core::WeightMode_CoordinateGaussian)
        {
            weight_token = 'C';
        }
        else if(weighted_options_.mode == weighted_moment_core::WeightMode_ValueGaussian)
        {
            weight_token = 'V';
        }
        else
        {
            weight_token = 'W';
        }
    }
    const int output_dim = (tempGridShape != NULL && tempGridShape->nx > 1 && tempGridShape->ny > 1 && tempGridShape->nz > 1) ? 3 : 2;
    const int order = Key != NULL ? Key->ORDER : 0;

    std::ostringstream name;
    name << (window_role ? "myWM" : "myGrid")
         << target_token
         << order
         << algorithm_token
         << weight_token
         << output_dim;
    return sanitize_grid_name_base(name.str());
}









int HOSC_TEMPLATE::set_MyGrid_Data(double** targetGrid, int gridSize, int valuePos)
{
	hosc_debug_log("set_MyGrid_Data.begin", "targetGrid=%p gridSize=%d valuePos=%d auto_grid=%p auto_property=%p", targetGrid, gridSize, valuePos, auto_grid_, auto_property_);
	if(auto_property_ == NULL || auto_grid_ == NULL || targetGrid == NULL)
	{
		GsTLcerr << "Cannot write HOSC output: target grid, target property, or result matrix is NULL." << gstlIO::end;
		hosc_debug_log("set_MyGrid_Data.fail", "null input");
		return -1;
	}
	if(valuePos < 0)
	{
		GsTLcerr << "Cannot write HOSC output: invalid value column." << gstlIO::end;
		hosc_debug_log("set_MyGrid_Data.fail", "invalid valuePos=%d", valuePos);
		return -1;
	}
	const int output_size = auto_grid_->nxyz();
	hosc_debug_log("set_MyGrid_Data.range", "rows=%d auto_nxyz=%d", gridSize, output_size);
	for(int node = 0; node < output_size; ++node)
	{
		auto_property_->set_not_informed(node);
	}
	std::vector<unsigned char> assigned((size_t)output_size, 0);
	for(int i = 0; i < gridSize; i++)
	{
		if(targetGrid[i] == NULL)
		{
			GsTLcerr << "Cannot write HOSC output: result row is NULL." << gstlIO::end;
			hosc_debug_log("set_MyGrid_Data.fail", "targetGrid[%d] is NULL", i);
			return -1;
		}
		const int node_id = template_result_row_to_node_id(output_layout_, auto_grid_, targetGrid[i]);
		if(node_id < 0 || node_id >= output_size || assigned[(size_t)node_id] != 0)
		{
			GsTLcerr << "Cannot write HOSC output: result lag coordinates do not map uniquely to the output grid." << gstlIO::end;
			hosc_debug_log("set_MyGrid_Data.fail", "row=%d node=%d lag mapping invalid", i, node_id);
			return -1;
		}
		assigned[(size_t)node_id] = 1;
		if(i == 0 || i == gridSize - 1)
		{
			hosc_debug_log("set_MyGrid_Data.sample", "row=%d node=%d value=%.17g", i, node_id, targetGrid[i][valuePos]);
		}
		const double value = targetGrid[i][valuePos];
		if(is_template_float_representable(value) && value != myEmpty)
		{
			auto_property_->set_value(value, node_id);
		}

	}


	hosc_debug_log("set_MyGrid_Data.end", "done");
	return 0;
}

bool HOSC_TEMPLATE::write_Template_Output_Validity(double** targetGrid, int gridSize, int valuePos)
{
	if(auto_grid_ == NULL || targetGrid == NULL || valuePos < 0)
	{
		return false;
	}
	const std::string numeric_property_name = auto_prop_name_ + "_numeric_valid";
	const std::string display_property_name = auto_prop_name_ + "_display_valid";
	GsTLGridProperty* numeric_validity = geostat_utils::add_property_to_grid(auto_grid_, numeric_property_name);
	GsTLGridProperty* display_validity = geostat_utils::add_property_to_grid(auto_grid_, display_property_name);
	if(numeric_validity == NULL || display_validity == NULL)
	{
		hosc_debug_log("write_Template_Output_Validity.fail", "property creation failed");
		return false;
	}
	const int output_size = auto_grid_->nxyz();
	for(int i = 0; i < output_size; ++i)
	{
		numeric_validity->set_value(0.0, i);
		display_validity->set_value(0.0, i);
	}
	for(int i = 0; i < gridSize; ++i)
	{
		if(targetGrid[i] == NULL)
		{
			return false;
		}
		const int node_id = template_result_row_to_node_id(output_layout_, auto_grid_, targetGrid[i]);
		if(node_id < 0 || node_id >= output_size)
		{
			return false;
		}
		const double value = targetGrid[i][valuePos];
		const bool numeric_valid = weighted_moment_core::is_finite_number(value) && value != myEmpty;
		numeric_validity->set_value(numeric_valid ? 1.0 : 0.0, node_id);
		display_validity->set_value(numeric_valid && is_template_float_representable(value) ? 1.0 : 0.0, node_id);
	}
	return true;
}

bool HOSC_TEMPLATE::write_Template_Exact_Support(double** targetGrid, int gridSize, int valuePos)
{
	if(auto_grid_ == NULL || targetGrid == NULL || valuePos < 0)
	{
		return false;
	}
	GsTLGridProperty* support = geostat_utils::add_property_to_grid(auto_grid_, auto_prop_name_ + "_support");
	if(support == NULL)
	{
		return false;
	}
	const int output_size = auto_grid_->nxyz();
	for(int i = 0; i < output_size; ++i)
	{
		support->set_value(0.0, i);
	}
	for(int i = 0; i < gridSize; ++i)
	{
		if(targetGrid[i] == NULL)
		{
			return false;
		}
		const int node_id = template_result_row_to_node_id(output_layout_, auto_grid_, targetGrid[i]);
		if(node_id < 0 || node_id >= output_size)
		{
			return false;
		}
		support->set_value(targetGrid[i][valuePos + 1], node_id);
	}
	return true;
}

RGrid* HOSC_TEMPLATE::CreateCartesianGridWithProperty(
	GridShape* tempGridShape,
	const char* gridName,
	const std::string& propertyName,
	GsTLGridProperty** out_property,
	std::string* created_name)
{
	if(out_property) *out_property = NULL;
	if(created_name) created_name->clear();
	hosc_debug_log(
		"CreateCartesianGridWithProperty.begin",
		"gridName=%s shape=%d,%d,%d property=%s",
		gridName == NULL ? "(null)" : gridName,
		tempGridShape != NULL ? tempGridShape->nx : -1,
		tempGridShape != NULL ? tempGridShape->ny : -1,
		tempGridShape != NULL ? tempGridShape->nz : -1,
		propertyName.c_str());
	if(current_project_ == NULL || tempGridShape == NULL)
	{
		hosc_debug_log("CreateCartesianGridWithProperty.fail", "invalid input");
		return NULL;
	}

	std::string base_name = gridName == NULL ? "SliceGrid" : gridName;
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
		hosc_debug_log("CreateCartesianGridWithProperty.fail", "failed to create unique name");
		return NULL;
	}

	SmartPtr<Named_interface> ni = Root::instance()->new_interface("cgrid", full_path);
	if(ni.raw_ptr() == NULL)
	{
		hosc_debug_log("CreateCartesianGridWithProperty.fail", "Root::new_interface failed");
		return NULL;
	}
	Cartesian_grid* grid = dynamic_cast<Cartesian_grid*>(ni.raw_ptr());
	if(grid == NULL)
	{
		hosc_debug_log("CreateCartesianGridWithProperty.fail", "dynamic_cast<Cartesian_grid*> failed");
		return NULL;
	}

	const int safe_nx = tempGridShape->nx > 0 ? tempGridShape->nx : 1;
	const int safe_ny = tempGridShape->ny > 0 ? tempGridShape->ny : 1;
	const int safe_nz = tempGridShape->nz > 0 ? tempGridShape->nz : 1;
	const double safe_dx = tempGridShape->dx > 0 ? tempGridShape->dx : 1.0;
	const double safe_dy = tempGridShape->dy > 0 ? tempGridShape->dy : 1.0;
	const double safe_dz = tempGridShape->dz > 0 ? tempGridShape->dz : 1.0;
	grid->set_dimensions(safe_nx, safe_ny, safe_nz, safe_dx, safe_dy, safe_dz);
	grid->origin(GsTLPoint(tempGridShape->cx, tempGridShape->cy, tempGridShape->cz));

	unsigned long register_exception = 0;
	hosc_debug_log("CreateCartesianGridWithProperty.register.begin", "unique=%s length=%d", unique_name.c_str(), (int)unique_name.size());
	if(!safe_register_project_object(current_project_, unique_name, &register_exception))
	{
		hosc_debug_log("CreateCartesianGridWithProperty.fail", "safe_register_project_object failed exception=0x%08lx", register_exception);
		return NULL;
	}

	SmartPtr<Named_interface> managed_ni = Root::instance()->interface(full_path);
	if(!managed_ni.raw_ptr())
	{
		hosc_debug_log("CreateCartesianGridWithProperty.fail", "managed grid not found");
		return NULL;
	}
	RGrid* managed_grid = dynamic_cast<RGrid*>(managed_ni.raw_ptr());
	Cartesian_grid* managed_cgrid = dynamic_cast<Cartesian_grid*>(managed_ni.raw_ptr());
	if(managed_grid == NULL || managed_cgrid == NULL)
	{
		hosc_debug_log("CreateCartesianGridWithProperty.fail", "managed cast failed");
		return NULL;
	}
	auto_property_ = geostat_utils::add_property_to_grid(managed_grid, propertyName);
	if(auto_property_ == NULL)
	{
		hosc_debug_log("CreateCartesianGridWithProperty.fail", "add_property_to_grid failed");
		return NULL;
	}
	managed_grid->set_level(1);
	if(out_property) *out_property = auto_property_;
	if(created_name) *created_name = unique_name;
	hosc_debug_log(
		"CreateCartesianGridWithProperty.end",
		"grid=%p property=%p name=%s size=%d",
		managed_grid,
		auto_property_,
		unique_name.c_str(),
		managed_grid->nxyz());
	return managed_grid;
}

bool HOSC_TEMPLATE::create_3d_slice_output()
{
	if(!slice_output_enabled_)
	{
		hosc_debug_log("create_3d_slice_output.skip", "disabled");
		return true;
	}
	if(auto_grid_ == NULL || auto_property_ == NULL)
	{
		hosc_debug_log("create_3d_slice_output.fail", "auto grid/property null");
		return false;
	}
	if(auto_grid_->nz() <= 1 || auto_grid_->nx() <= 1 || auto_grid_->ny() <= 1)
	{
		hosc_debug_log(
			"create_3d_slice_output.skip",
			"not a 3d cartesian output nx=%d ny=%d nz=%d",
			auto_grid_->nx(),
			auto_grid_->ny(),
			auto_grid_->nz());
		return true;
	}
	slice_grid_ = NULL;
	slice_property_ = NULL;
	slice_grid_name_.clear();

	RGrid* source_grid = auto_grid_;
	GsTLGridProperty* source_property = auto_property_;
	const int plane = slice_plane_mode_;
	const int pos = slice_position_;
	int slice_index = pos - 1;
	int slice_nx = 1;
	int slice_ny = 1;
	int slice_nz = 1;
	GridShape slice_shape = myGridShape;
	std::string slice_name = auto_grid_name_.empty() ? "myGrid" : auto_grid_name_;
	std::ostringstream suffix;

	const int max_x = source_grid->nx();
	const int max_y = source_grid->ny();
	const int max_z = source_grid->nz();

	double source_min = 0.0;
	double source_max = 0.0;
	double source_abs_sum = 0.0;
	int source_nonzero = 0;
	const int source_count = source_grid->nxyz();
	for(int i = 0; i < source_count; i++)
	{
		const double value = source_property->get_value(i);
		if(i == 0 || value < source_min)
		{
			source_min = value;
		}
		if(i == 0 || value > source_max)
		{
			source_max = value;
		}
		if(value != 0.0)
		{
			source_nonzero++;
		}
		source_abs_sum += fabs(value);
	}
	hosc_debug_log(
		"create_3d_slice_output.source_stats",
		"grid=%s property=%s shape=%d,%d,%d count=%d nonzero=%d min=%.17g max=%.17g sum_abs=%.17g",
		slice_name.c_str(),
		auto_prop_name_.c_str(),
		max_x,
		max_y,
		max_z,
		source_count,
		source_nonzero,
		source_min,
		source_max,
		source_abs_sum);
	for(int z = 0; z < max_z; z++)
	{
		double layer_min = 0.0;
		double layer_max = 0.0;
		double layer_abs_sum = 0.0;
		int layer_nonzero = 0;
		const int layer_count = max_x * max_y;
		for(int y = 0; y < max_y; y++)
		{
			for(int x = 0; x < max_x; x++)
			{
				const int src_id = x + y * max_x + z * max_x * max_y;
				const double value = source_property->get_value(src_id);
				if((x == 0 && y == 0) || value < layer_min)
				{
					layer_min = value;
				}
				if((x == 0 && y == 0) || value > layer_max)
				{
					layer_max = value;
				}
				if(value != 0.0)
				{
					layer_nonzero++;
				}
				layer_abs_sum += fabs(value);
			}
		}
		hosc_debug_log(
			"create_3d_slice_output.xy_layer_stats",
			"z_pos=%d count=%d nonzero=%d min=%.17g max=%.17g sum_abs=%.17g",
			z + 1,
			layer_count,
			layer_nonzero,
			layer_min,
			layer_max,
			layer_abs_sum);
	}
	if(plane == SlicePlane_XY)
	{
		if(pos < 1 || pos > max_z)
		{
			hosc_debug_log("create_3d_slice_output.fail", "xy position out of range pos=%d max_z=%d", pos, max_z);
			return false;
		}
		slice_nx = max_x;
		slice_ny = max_y;
		slice_nz = 1;
		slice_shape.nx = slice_nx;
		slice_shape.ny = slice_ny;
		slice_shape.nz = 1;
		slice_shape.cz = source_grid->location(0).z();
		suffix << "_xy_z" << pos;
	}
	else if(plane == SlicePlane_XZ)
	{
		if(pos < 1 || pos > max_y)
		{
			hosc_debug_log("create_3d_slice_output.fail", "xz position out of range pos=%d max_y=%d", pos, max_y);
			return false;
		}
		slice_nx = max_x;
		slice_ny = max_z;
		slice_nz = 1;
		slice_shape.nx = slice_nx;
		slice_shape.ny = slice_ny;
		slice_shape.nz = 1;
		suffix << "_xz_y" << pos;
	}
	else if(plane == SlicePlane_YZ)
	{
		if(pos < 1 || pos > max_x)
		{
			hosc_debug_log("create_3d_slice_output.fail", "yz position out of range pos=%d max_x=%d", pos, max_x);
			return false;
		}
		slice_nx = max_y;
		slice_ny = max_z;
		slice_nz = 1;
		slice_shape.nx = slice_nx;
		slice_shape.ny = slice_ny;
		slice_shape.nz = 1;
		suffix << "_yz_x" << pos;
	}
	else
	{
		hosc_debug_log("create_3d_slice_output.fail", "invalid plane mode=%d", plane);
		return false;
	}

	std::string slice_prop_name = auto_prop_name_;
	if(slice_prop_name.empty())
	{
		slice_prop_name = "value";
	}
	slice_prop_name += suffix.str();
	GsTLGridProperty* slice_property = NULL;
	std::string created_slice_name;
	RGrid* slice_grid = CreateCartesianGridWithProperty(
		&slice_shape,
		(slice_name + suffix.str()).c_str(),
		slice_prop_name,
		&slice_property,
		&created_slice_name);
	if(slice_grid == NULL || slice_property == NULL)
	{
		hosc_debug_log("create_3d_slice_output.fail", "slice grid create failed");
		return false;
	}

	const int count = slice_grid->nxyz();
	double slice_min = 0.0;
	double slice_max = 0.0;
	double slice_abs_sum = 0.0;
	int slice_nonzero = 0;
	for(int i = 0; i < count; i++)
	{
		const int x = i % slice_nx;
		const int y = (i / slice_nx) % slice_ny;
		int src_id = -1;
		if(plane == SlicePlane_XY)
		{
			src_id = slice_index * max_x * max_y + y * max_x + x;
		}
		else if(plane == SlicePlane_XZ)
		{
			const int z = y;
			src_id = x + slice_index * max_x + z * max_x * max_y;
		}
		else
		{
			const int z = y;
			const int source_y = x;
			src_id = slice_index + source_y * max_x + z * max_x * max_y;
		}
		if(src_id < 0 || src_id >= source_grid->nxyz())
		{
			hosc_debug_log("create_3d_slice_output.fail", "invalid source node id at i=%d src=%d", i, src_id);
			return false;
		}
		const double value = source_property->get_value(src_id);
		if(i == 0 || value < slice_min)
		{
			slice_min = value;
		}
		if(i == 0 || value > slice_max)
		{
			slice_max = value;
		}
		if(value != 0.0)
		{
			slice_nonzero++;
		}
		slice_abs_sum += fabs(value);
		if(i < 5 || i == count - 1)
		{
			hosc_debug_log(
				"create_3d_slice_output.copy_sample",
				"dst=%d src=%d value=%.17g",
				i,
				src_id,
				value);
		}
		slice_property->set_value(value, i);
	}

	slice_grid_ = slice_grid;
	slice_property_ = slice_property;
	slice_grid_name_ = created_slice_name;
	hosc_debug_log(
		"create_3d_slice_output.slice_stats",
		"plane=%s position=%d count=%d nonzero=%d min=%.17g max=%.17g sum_abs=%.17g",
		slice_plane_label(plane),
		pos,
		count,
		slice_nonzero,
		slice_min,
		slice_max,
		slice_abs_sum);
	hosc_debug_log(
		"create_3d_slice_output.end",
		"plane=%s position=%d output_name=%s size=%d",
		slice_plane_label(plane),
		pos,
		slice_grid_name_.c_str(),
		count);
	return true;
}

bool HOSC_TEMPLATE::create_3d_slice_output_from_selected_grid()
{
	if(!slice_output_enabled_)
	{
		hosc_debug_log("create_3d_slice_output_from_selected_grid.skip", "disabled");
		return true;
	}
	if(current_project_ == NULL)
	{
		hosc_debug_log("create_3d_slice_output_from_selected_grid.fail", "current_project_ null");
		return false;
	}
	if(slice_source_grid_name_.empty())
	{
		GsTLcerr << "No 3D slice source grid selected." << gstlIO::end;
		hosc_debug_log("create_3d_slice_output_from_selected_grid.fail", "empty source grid");
		return false;
	}

	SmartPtr<Named_interface> source_ni =
		Root::instance()->interface(gridModels_manager + "/" + slice_source_grid_name_);
	RGrid* source_rgrid = dynamic_cast<RGrid*>(source_ni.raw_ptr());
	Geostat_grid* source_geostat_grid = dynamic_cast<Geostat_grid*>(source_ni.raw_ptr());
	if(source_rgrid == NULL || source_geostat_grid == NULL)
	{
		GsTLcerr << "Selected 3D slice source is not a Cartesian grid." << gstlIO::end;
		hosc_debug_log("create_3d_slice_output_from_selected_grid.fail", "source lookup/cast failed grid=%s", slice_source_grid_name_.c_str());
		return false;
	}
	if(source_rgrid->nx() <= 1 || source_rgrid->ny() <= 1 || source_rgrid->nz() <= 1)
	{
		GsTLcerr << "Selected slice source must be a three-dimensional Cartesian grid." << gstlIO::end;
		hosc_debug_log(
			"create_3d_slice_output_from_selected_grid.fail",
			"not 3d nx=%d ny=%d nz=%d",
			source_rgrid->nx(),
			source_rgrid->ny(),
			source_rgrid->nz());
		return false;
	}

	GsTLGridProperty* source_property = NULL;
	if(!slice_source_property_name_.empty())
	{
		source_property = source_geostat_grid->property(slice_source_property_name_);
	}
	if(source_property == NULL)
	{
		source_property = source_geostat_grid->selected_property();
	}
	if(source_property == NULL)
	{
		std::list<std::string> properties = source_geostat_grid->property_list();
		if(!properties.empty())
		{
			source_property = source_geostat_grid->property(properties.front());
		}
	}
	if(source_property == NULL)
	{
		GsTLcerr << "Selected 3D slice source has no readable property." << gstlIO::end;
		hosc_debug_log("create_3d_slice_output_from_selected_grid.fail", "source property null");
		return false;
	}

	auto_grid_ = source_rgrid;
	auto_property_ = source_property;
	auto_grid_name_ = slice_source_grid_name_;
	auto_prop_name_ = source_property->name();
	myGridShape.nx = source_rgrid->nx();
	myGridShape.ny = source_rgrid->ny();
	myGridShape.nz = source_rgrid->nz();
	myGridShape.dx = 1.0;
	myGridShape.dy = 1.0;
	myGridShape.dz = 1.0;
	if(source_rgrid->nxyz() > 0)
	{
		GsTLPoint p0 = source_rgrid->location(0);
		myGridShape.cx = p0.x();
		myGridShape.cy = p0.y();
		myGridShape.cz = p0.z();
	}
	else
	{
		myGridShape.cx = 0.0;
		myGridShape.cy = 0.0;
		myGridShape.cz = 0.0;
	}

	hosc_debug_log(
		"create_3d_slice_output_from_selected_grid.begin",
		"source=%s property=%s shape=%d,%d,%d plane=%s pos=%d",
		auto_grid_name_.c_str(),
		auto_prop_name_.c_str(),
		myGridShape.nx,
		myGridShape.ny,
		myGridShape.nz,
		slice_plane_label(slice_plane_mode_),
		slice_position_);
	return create_3d_slice_output();
}

bool HOSC_TEMPLATE::get_My_Shape(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	output_layout_.valid = false;
	if(Key == NULL)
	{
		error_mesgs->report("template_shape", "Keys are not initialized");
		hosc_debug_log("get_My_Shape.fail", "Key is NULL");
		return false;
	}

	hosc_debug_log("get_My_Shape.begin", "start generalized cross-section layout");
	auto_prop_name_ = parameters->value("props.value");
	if(auto_prop_name_.empty())
	{
		error_mesgs->report("props", "No property specified");
		hosc_debug_log("get_My_Shape.fail", "empty output property");
		return false;
	}

	const double dcx = atof(parameters->value("xmn_.value").c_str());
	const double dcy = atof(parameters->value("ymn_.value").c_str());
	const double dcz = atof(parameters->value("zmn_.value").c_str());

	const int raw_dirs[3] = {
		atoi(parameters->value("fst_dir.value").c_str()),
		atoi(parameters->value("snd_dir.value").c_str()),
		atoi(parameters->value("trd_dir.value").c_str())
	};
	std::vector<int> displayedDirs;
	for(int i = 0; i < 3; i++)
	{
		const int dir = raw_dirs[i];
		if(dir == 0)
		{
			continue;
		}
		if(dir < 1 || dir > Key->NDIR)
		{
			error_mesgs->report(i == 0 ? "fst_dir" : (i == 1 ? "snd_dir" : "trd_dir"), "Direction index is out of range");
			hosc_debug_log("get_My_Shape.fail", "dir[%d]=%d out of range ndir=%d", i, dir, Key->NDIR);
			return false;
		}
		for(size_t j = 0; j < displayedDirs.size(); j++)
		{
			if(displayedDirs[j] == dir)
			{
				error_mesgs->report(i == 0 ? "fst_dir" : (i == 1 ? "snd_dir" : "trd_dir"), "Displayed directions must be distinct");
				hosc_debug_log("get_My_Shape.fail", "duplicate displayed direction=%d", dir);
				return false;
			}
		}
		displayedDirs.push_back(dir);
	}

	if(displayedDirs.size() != 2 && displayedDirs.size() != 3)
	{
		error_mesgs->report("fst_dir", "Select exactly two directions for 2D output or three directions for 3D output");
		hosc_debug_log("get_My_Shape.fail", "displayed_count=%d", (int)displayedDirs.size());
		return false;
	}

	std::vector<int> hiddenDirs;
	for(int dir = 1; dir <= Key->NDIR; dir++)
	{
		bool isDisplayed = false;
		for(size_t i = 0; i < displayedDirs.size(); i++)
		{
			if(displayedDirs[i] == dir)
			{
				isDisplayed = true;
				break;
			}
		}
		if(!isDisplayed)
		{
			hiddenDirs.push_back(dir);
		}
	}

	std::vector<int> fixedPositions;
	const std::string crossectionContent = parameters->value("cross_pos.value");
	hosc_debug_log("get_My_Shape.cross_raw", "cross_pos=%s hidden_count=%d", crossectionContent.c_str(), (int)hiddenDirs.size());
	if(crossectionContent.empty())
	{
		fixedPositions.assign(hiddenDirs.size(), 1);
	}
	else
	{
		if(crossectionContent[0] == ',' || crossectionContent[crossectionContent.size() - 1] == ',')
		{
			error_mesgs->report("cross_pos", "Cross-section positions must not begin or end with comma");
			hosc_debug_log("get_My_Shape.fail", "cross_pos begins/ends with comma");
			return false;
		}
		for(size_t i = 0; i < crossectionContent.size(); i++)
		{
			const char ch = crossectionContent[i];
			if(!(isdigit((unsigned char)ch) || ch == ','))
			{
				error_mesgs->report("cross_pos", "Cross-section positions must be positive integers separated by commas");
				hosc_debug_log("get_My_Shape.fail", "invalid cross_pos char index=%d char=%d", (int)i, (int)ch);
				return false;
			}
			if(ch == ',' && i + 1 < crossectionContent.size() && crossectionContent[i + 1] == ',')
			{
				error_mesgs->report("cross_pos", "Cross-section positions contain an empty item");
				hosc_debug_log("get_My_Shape.fail", "double comma in cross_pos");
				return false;
			}
		}

		std::stringstream ss(crossectionContent);
		std::string item;
		while(std::getline(ss, item, ','))
		{
			if(item.empty())
			{
				error_mesgs->report("cross_pos", "Cross-section positions contain an empty item");
				hosc_debug_log("get_My_Shape.fail", "empty cross_pos item");
				return false;
			}
			const int pos = atoi(item.c_str());
			if(pos <= 0)
			{
				error_mesgs->report("cross_pos", "Cross-section positions must be positive");
				hosc_debug_log("get_My_Shape.fail", "non-positive cross value=%d", pos);
				return false;
			}
			fixedPositions.push_back(pos);
		}
	}

	if(hiddenDirs.empty())
	{
		if(fixedPositions.size() > 1 || (fixedPositions.size() == 1 && fixedPositions[0] != 1))
		{
			error_mesgs->report("cross_pos", "No fixed cross-section position is needed when all directions are displayed");
			hosc_debug_log("get_My_Shape.fail", "unexpected fixed positions count=%d", (int)fixedPositions.size());
			return false;
		}
		fixedPositions.clear();
	}
	else if(fixedPositions.size() != hiddenDirs.size())
	{
		error_mesgs->report("cross_pos", "The number of positions must equal the number of hidden directions");
		hosc_debug_log("get_My_Shape.fail", "cross_count=%d expected=%d", (int)fixedPositions.size(), (int)hiddenDirs.size());
		return false;
	}

	for(size_t i = 0; i < fixedPositions.size(); i++)
	{
		const int originalDir = hiddenDirs[i];
		if(fixedPositions[i] < 1 || fixedPositions[i] > Key->NLAG[originalDir - 1])
		{
			error_mesgs->report("cross_pos", "Cross-section position is out of range for its hidden direction");
			hosc_debug_log(
				"get_My_Shape.fail",
				"cross[%d]=%d hidden_dir=%d original_nlag=%d",
				(int)i,
				fixedPositions[i],
				originalDir,
				Key->NLAG[originalDir - 1]);
			return false;
		}
	}

	myGridShape.nx = Key->NLAG[displayedDirs[0] - 1];
	myGridShape.ny = Key->NLAG[displayedDirs[1] - 1];
	myGridShape.nz = displayedDirs.size() == 3 ? Key->NLAG[displayedDirs[2] - 1] : 1;
	myGridShape.cx = dcx;
	myGridShape.cy = dcy;
	myGridShape.cz = dcz;
	myGridShape.dx = 1;
	myGridShape.dy = 1;
	myGridShape.dz = 1;

	std::vector<int> oldNOD(Key->NOD, Key->NOD + Key->NDIR);
	std::vector<int> oldNLAG(Key->NLAG, Key->NLAG + Key->NDIR);
	std::vector<int> oldLAG(Key->LAG, Key->LAG + Key->NDIR);
	std::vector<int> oldLTOL(Key->LTOL, Key->LTOL + Key->NDIR);
	std::vector<int> oldAZM(Key->AZM, Key->AZM + Key->NDIR);
	std::vector<int> oldATOL(Key->ATOL, Key->ATOL + Key->NDIR);
	std::vector<int> oldBANDWH(Key->BANDWH, Key->BANDWH + Key->NDIR);
	std::vector<int> oldDIP(Key->DIP, Key->DIP + Key->NDIR);
	std::vector<int> oldDTOL(Key->DTOL, Key->DTOL + Key->NDIR);
	std::vector<int> oldBANDWD(Key->BANDWD, Key->BANDWD + Key->NDIR);

	std::vector<int> reorderedDirs = hiddenDirs;
	reorderedDirs.insert(reorderedDirs.end(), displayedDirs.begin(), displayedDirs.end());
	std::vector<ExactTemplateDirection> oldExactDirections = integer_direction_vectors_;
	for(int i = 0; i < Key->NDIR; i++)
	{
		const int oldIndex = reorderedDirs[i] - 1;
		Key->NOD[i] = oldNOD[oldIndex];
		Key->NLAG[i] = oldNLAG[oldIndex];
		Key->LAG[i] = oldLAG[oldIndex];
		Key->LTOL[i] = oldLTOL[oldIndex];
		Key->AZM[i] = oldAZM[oldIndex];
		Key->ATOL[i] = oldATOL[oldIndex];
		Key->BANDWH[i] = oldBANDWH[oldIndex];
		Key->DIP[i] = oldDIP[oldIndex];
		Key->DTOL[i] = oldDTOL[oldIndex];
		Key->BANDWD[i] = oldBANDWD[oldIndex];
		if(template_sampling_exact_ && oldIndex >= 0 && oldIndex < (int)oldExactDirections.size())
		{
			integer_direction_vectors_[i] = oldExactDirections[oldIndex];
		}
	}

	crossectionNum = (int)fixedPositions.size();
	for(int i = 0; i < 10; i++)
	{
		crossection[i] = 0;
	}
	if(crossectionNum > 10)
	{
		error_mesgs->report("cross_pos", "Too many fixed cross-section positions");
		hosc_debug_log("get_My_Shape.fail", "cross_count=%d exceeds storage", crossectionNum);
		return false;
	}
	for(int i = 0; i < crossectionNum; i++)
	{
		crossection[i] = fixedPositions[i] - 1;
	}

	// The calculators call Destruction() before returning targetGrid. Freeze the
	// reordered result-row columns and lag steps before that ownership boundary.
	output_layout_.displayed_direction_count = (int)displayedDirs.size();
	output_layout_.dimensions[0] = myGridShape.nx;
	output_layout_.dimensions[1] = myGridShape.ny;
	output_layout_.dimensions[2] = myGridShape.nz;
	for(int axis = 0; axis < 3; ++axis)
	{
		output_layout_.row_column[axis] = -1;
		output_layout_.lag_step[axis] = 0;
	}
	for(int axis = 0; axis < output_layout_.displayed_direction_count; ++axis)
	{
		const int reordered_direction = Key->NDIR - output_layout_.displayed_direction_count + axis;
		if(reordered_direction < 0 || reordered_direction >= Key->NDIR || Key->LAG == NULL || Key->LAG[reordered_direction] <= 0)
		{
			error_mesgs->report("template_shape", "Unable to freeze the displayed lag layout");
			hosc_debug_log("get_My_Shape.fail", "cannot freeze layout axis=%d reordered_direction=%d", axis, reordered_direction);
			return false;
		}
		output_layout_.row_column[axis] = reordered_direction;
		output_layout_.lag_step[axis] = Key->LAG[reordered_direction];
	}
	if(output_layout_.dimensions[0] <= 0 || output_layout_.dimensions[1] <= 0 ||
	   output_layout_.dimensions[2] <= 0 ||
	   (output_layout_.displayed_direction_count == 2 && output_layout_.dimensions[2] != 1))
	{
		error_mesgs->report("template_shape", "The frozen template output layout is invalid");
		hosc_debug_log("get_My_Shape.fail", "invalid frozen layout count=%d dimensions=%d,%d,%d", output_layout_.displayed_direction_count, output_layout_.dimensions[0], output_layout_.dimensions[1], output_layout_.dimensions[2]);
		return false;
	}
	output_layout_.valid = true;

	hosc_debug_log(
		"get_My_Shape.end",
		"shape=%d,%d,%d displayed=%d,%d,%d hidden_count=%d cross_count=%d cross0_zero=%d layout_count=%d layout_columns=%d,%d,%d layout_steps=%d,%d,%d",
		myGridShape.nx,
		myGridShape.ny,
		myGridShape.nz,
		displayedDirs.size() > 0 ? displayedDirs[0] : 0,
		displayedDirs.size() > 1 ? displayedDirs[1] : 0,
		displayedDirs.size() > 2 ? displayedDirs[2] : 0,
		(int)hiddenDirs.size(),
		crossectionNum,
		crossectionNum > 0 ? crossection[0] : -1,
		output_layout_.displayed_direction_count,
		output_layout_.row_column[0],
		output_layout_.row_column[1],
		output_layout_.row_column[2],
		output_layout_.lag_step[0],
		output_layout_.lag_step[1],
		output_layout_.lag_step[2]);
	return true;

}

bool HOSC_TEMPLATE::get_Output_Slice_Options(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	slice_output_enabled_ = parse_bool_selection(parameter_value_or_empty(parameters, "Output_3D_Slice_Enable.value"));
	slice_plane_mode_ = parse_slice_plane_selection(parameter_value_or_empty(parameters, "Output_3D_Slice_Plane.value"));
	slice_position_ = atoi(parameter_value_or_empty(parameters, "Output_3D_Slice_Position.value").c_str());

	if(!slice_output_enabled_)
	{
		hosc_debug_log("get_Output_Slice_Options.end", "disabled");
		return true;
	}

	if(slice_plane_mode_ < 0)
	{
		error_mesgs->report("Output_3D_Slice_Plane", "Slice plane must be xy, xz, or yz");
		hosc_debug_log("get_Output_Slice_Options.fail", "invalid plane");
		return false;
	}
	if(slice_position_ <= 0)
	{
		error_mesgs->report("Output_3D_Slice_Position", "Slice position must be positive");
		hosc_debug_log("get_Output_Slice_Options.fail", "invalid position=%d", slice_position_);
		return false;
	}
	hosc_debug_log(
		"get_Output_Slice_Options.end",
		"enabled=%d plane=%s position=%d",
		slice_output_enabled_ ? 1 : 0,
		slice_plane_label(slice_plane_mode_),
		slice_position_);
	return true;
}

bool HOSC_TEMPLATE::get_Output_Slice_Target(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	if(!slice_output_enabled_)
	{
		return true;
	}
	slice_source_grid_name_ = parameter_value_or_empty(parameters, "Output_3D_Slice_Source_Grid.value");
	slice_source_property_name_ = parameter_value_or_empty(parameters, "Output_3D_Slice_Source_Property.value");

	if(slice_source_grid_name_.empty() || slice_source_grid_name_ == "<- None ->")
	{
		error_mesgs->report("Output_3D_Slice_Source_Grid", "Select a 3D Cartesian grid to slice");
		hosc_debug_log("get_Output_Slice_Target.fail", "empty source grid");
		return false;
	}

	SmartPtr<Named_interface> source_ni =
		Root::instance()->interface(gridModels_manager + "/" + slice_source_grid_name_);
	RGrid* source_rgrid = dynamic_cast<RGrid*>(source_ni.raw_ptr());
	Geostat_grid* source_geostat_grid = dynamic_cast<Geostat_grid*>(source_ni.raw_ptr());
	if(source_rgrid == NULL || source_geostat_grid == NULL)
	{
		error_mesgs->report("Output_3D_Slice_Source_Grid", "The selected slice source must be a Cartesian grid");
		hosc_debug_log("get_Output_Slice_Target.fail", "source cast failed grid=%s", slice_source_grid_name_.c_str());
		return false;
	}
	if(source_rgrid->nx() <= 1 || source_rgrid->ny() <= 1 || source_rgrid->nz() <= 1)
	{
		error_mesgs->report("Output_3D_Slice_Source_Grid", "The selected slice source must be three-dimensional");
		hosc_debug_log(
			"get_Output_Slice_Target.fail",
			"not 3d nx=%d ny=%d nz=%d",
			source_rgrid->nx(),
			source_rgrid->ny(),
			source_rgrid->nz());
		return false;
	}
	if(!slice_source_property_name_.empty() && source_geostat_grid->property(slice_source_property_name_) == NULL)
	{
		error_mesgs->report("Output_3D_Slice_Source_Property", "The selected source property does not exist");
		hosc_debug_log("get_Output_Slice_Target.fail", "property not found=%s", slice_source_property_name_.c_str());
		return false;
	}
	if(slice_source_property_name_.empty() && source_geostat_grid->selected_property() == NULL && source_geostat_grid->property_list().empty())
	{
		error_mesgs->report("Output_3D_Slice_Source_Property", "The selected source grid has no property");
		hosc_debug_log("get_Output_Slice_Target.fail", "no property");
		return false;
	}

	hosc_debug_log(
		"get_Output_Slice_Target.end",
		"grid=%s property=%s",
		slice_source_grid_name_.c_str(),
		slice_source_property_name_.c_str());
	return true;
}

#if 0
RGrid* HOSC_TEMPLATE::getDefaultGrid()
{
    std::string grid_path = gridModels_manager + "/HOSC_DefaultGrid";
    SmartPtr<Named_interface> ni = Root::instance()->interface(grid_path);
    return dynamic_cast<RGrid*>(ni.raw_ptr());
}

bool HOSC_TEMPLATE::get_Default_Data(const Parameters_handler* parameters, Error_messages_handler* error_mesgs)
{
	if(default_grid_ == NULL)
		return false;
	string default_prop_name_ = parameters->value( "props.value" );
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
            GsTLlog << "HOSC_TEMPLATE Plugin initializing..." << gstlIO::end;


            SmartPtr<Named_interface> ui_ni =
                Root::instance()->new_interface(
                    "HOSC_TEMPLATE.ui",
                    "HOSC_TEMPLATE"
                );

            if (ui_ni.raw_ptr()) {
                GsTLlog << "HOSC_TEMPLATE UI registered successfully" << gstlIO::end;
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


GEOSTAT_PLUGIN(HOSC_TEMPLATE);
