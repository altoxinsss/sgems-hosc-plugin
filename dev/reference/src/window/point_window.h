#ifndef WITHOUT_SGEMS_POINT_WINDOW_H
#define WITHOUT_SGEMS_POINT_WINDOW_H

#include "window_geometry.h"
#include "point_set.h"
#include "text.h"
#include "weight_config.h"
#include "window_weights.h"
#include "WeightedMomentCore.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace without_sgems {
namespace point_window_detail {

inline bool boolean(const std::string& value) {
    const std::string v = to_lower(value);
    if (v == "true" || v == "1" || v == "yes" || v == "on") return true;
    if (v == "false" || v == "0" || v == "no" || v == "off") return false;
    throw std::runtime_error("Invalid Boolean value: " + value);
}

inline int positive_integer(const std::string& value, const std::string& name) {
    size_t end = 0;
    int result = 0;
    try { result = std::stoi(value, &end); }
    catch (...) { throw std::runtime_error(name + " must be a positive integer."); }
    if (end != value.size() || result <= 0)
        throw std::runtime_error(name + " must be a positive integer.");
    return result;
}

inline double finite_number(const std::string& text, const std::string& name) {
    const std::string value = trim(text);
    size_t end = 0;
    double result;
    try { result = std::stod(value, &end); }
    catch (...) { throw std::runtime_error(name + " must contain a finite numeric value."); }
    if (end != value.size() || !std::isfinite(result))
        throw std::runtime_error(name + " must contain a finite numeric value, without trailing text.");
    return result;
}

inline std::vector<double> finite_list(const std::string& text, const std::string& name) {
    const std::string value = trim(text);
    std::vector<double> result;
    std::istringstream stream(value);
    std::string token;
    if (value.find(',') != std::string::npos) {
        if (value.back() == ',') throw std::runtime_error(name + " has an empty numeric field.");
        while (std::getline(stream, token, ',')) result.push_back(finite_number(token, name));
    }
    else {
        while (stream >> token) result.push_back(finite_number(token, name));
    }
    return result;
}

struct Evaluator {
    bool bell = false;
    bool weighted = false;
};

inline Evaluator parse_evaluator(const std::string& text) {
    std::string name = to_lower(trim(text));
    for (char& c : name) if (c == '_' || c == '-') c = ' ';
    std::istringstream words(name);
    std::string word;
    bool avl = false, bell = false, weighted = false;
    while (words >> word) {
        if (word == "avl") avl = true;
        else if (word == "bell") bell = true;
        else if (word == "weighted") weighted = true;
        else if (word != "based" && word != "algorithm" && word != "partition" && word != "reference")
            throw std::runtime_error("Unknown point-set calculate_algorithm: " + text + ". Use Bell or AVL Based algorithm.");
    }
    if ((!avl && !bell && !weighted) || (avl && bell))
        throw std::runtime_error("Point-set calculate_algorithm must select Bell or AVL, not conflicting routes.");
    // 'Weighted' alone is the historical uniform recurrence alias. An explicit
    // AVL selection retains AVL evaluation even with the weighted qualifier.
    return { !avl, weighted };
}

typedef hosc_window_geometry::Shape Ellipsoid;

inline Ellipsoid parse_ellipsoid(const IniConfig& config, int dimensions) {
    for (const char* key : {"source_origin", "source_spacing", "target_origin"})
        if (config.has("general", key)) throw std::runtime_error(std::string("Point-set coordinates are physical; do not use ") + key);
    for (const char* key : {"window_size", "windowl", "windoww", "windowh", "move_step",
            "stepx", "stepy", "stepz", "require_full_window", "require_complete_window",
            "require_anchor_present", "output_at_window_center"})
        if (config.has("sliding_window", key)) throw std::runtime_error(std::string("Point-set windows use physical radii and grid_spacing, not ") + key);
    const char* section="sliding_window";
    std::vector<double> radii,angles={90,0,0};
    if(config.has(section,"neigh_ellipsoid")){
        if(config.has(section,"neighborhood_radii") || config.has(section,"neighborhood_angles"))
            throw std::runtime_error("Use neigh_ellipsoid or neighborhood_radii/neighborhood_angles, not both.");
        const auto v=finite_list(config.get(section,"neigh_ellipsoid"),"neigh_ellipsoid");
        if(v.size()!=6)throw std::runtime_error("neigh_ellipsoid requires six numbers.");
        radii.assign(v.begin(),v.begin()+3);angles.assign(v.begin()+3,v.end());
        if(config.get_or(section,"neighborhood_shape","ellipsoid")!="ellipsoid")
            throw std::runtime_error("neigh_ellipsoid is only valid with ellipsoid shape.");
    }else{
        radii=finite_list(config.get(section,"neighborhood_radii"),"neighborhood_radii");
        angles=finite_list(config.get_or(section,"neighborhood_angles","90,0,0"),"neighborhood_angles");
    }
    return Ellipsoid(config.get_or(section,"neighborhood_shape","ellipsoid"),radii,angles,dimensions);
}
inline hosc_window_geometry::CenterGrid center_grid(const IniConfig& config,
    const PointSetData& data,const Ellipsoid& shape) {
    const char* section="sliding_window";
    hosc_window_geometry::CenterGrid grid;
    const std::string mode=config.get_or(section,"center_mode","sample_points");
    const std::string bounds=config.get_or(section,"grid_bounds_mode","auto");
    auto list=[&](const char* key) {return finite_list(config.get_or(section,key,""),key);};
    const std::vector<double> empty;
    grid.configure(data,shape,mode,bounds,
        mode=="regular_grid"&&bounds=="manual"?list("grid_bounds_min"):empty,
        mode=="regular_grid"&&bounds=="manual"?list("grid_bounds_max"):empty,
        mode=="regular_grid"?list("grid_start"):empty,
        mode=="regular_grid"?list("grid_spacing"):empty);
    return grid;
}

inline double percentile(std::vector<double> values, double pct) {
    std::sort(values.begin(), values.end());
    const double position = pct * 0.01 * static_cast<double>(values.size() - 1);
    const size_t lo = static_cast<size_t>(std::floor(position));
    const size_t hi = static_cast<size_t>(std::ceil(position));
    return values[lo] + (position - lo) * (values[hi] - values[lo]);
}

} // namespace point_window_detail

inline int run_point_window(
    const IniConfig& config, const std::filesystem::path& config_path,
    const std::filesystem::path& source_path, const std::string& property,
    int dimensions, const std::string& target, const std::string& statistic,
    const std::string& algorithm,
    const std::function<double(const std::vector<double>&, int, bool, const std::vector<double>&)>& bell,
    const std::function<double(const std::vector<double>&, int, bool, const std::vector<double>&)>& avl) {
    using namespace point_window_detail;
    const auto options = parse_checked_weight_options(config, true);
    if (config.has("general", "coordinate_dimensions")) {
        const std::string raw = trim(config.get("general", "coordinate_dimensions"));
        if (raw != "0" && raw != "1" && raw != "2" && raw != "3")
            throw std::runtime_error("coordinate_dimensions must be 0, 1, 2 or 3.");
    }
    const std::string statistic_name = to_lower(trim(statistic));
    if (statistic_name != "moment" && statistic_name != "cumulant")
        throw std::runtime_error("Point-set calculate_target must be Moment or Cumulant.");
    const Evaluator evaluator = parse_evaluator(algorithm);
    const double missing = config.has("general", "missing_value")
        ? finite_number(config.get("general", "missing_value"), "missing_value")
        : std::numeric_limits<double>::quiet_NaN();
    const PointSetData data = load_point_set(source_path, property, dimensions, missing, external_weight_property(config));
    const WindowWeights source_weights(config, data);
    const Ellipsoid neighborhood = parse_ellipsoid(config, data.dimensions);
    const auto centers = center_grid(config, data, neighborhood);
    const int order = positive_integer(config.get_or("sliding_window", "moment_cumulant_order",
        config.get_or("sliding_window", "moment_order", "0")), "Moment/Cumulant order");
    if (order > 64) throw std::runtime_error("Point-set cumulant order must be at most 64.");
    const int minimum = positive_integer(config.get_or("sliding_window", "min_samples",
        config.get_or("sliding_window", "minimum_samples", "1")), "min_samples");
    const bool includes_center = boolean(config.get_or("sliding_window", "includes_center", "true"));
    if (centers.regular && !includes_center)
        throw std::runtime_error("includes_center=false applies only to sample_points mode; regular centers have no source record identity.");
    const bool center = boolean(config.get_or("sliding_window", "center_cumulant_input_for_order_ge_2", "false"));
    const bool invert = boolean(config.get_or("sliding_window", "invert_window_output", "false"));
    const bool clip = boolean(config.get_or("sliding_window", "enable_percentile_clipping",
        config.get_or("sliding_window", "percentile_clipping", "false")));
    double clip_low_pct = finite_number(config.get_or("sliding_window", "percentile_vmin",
        config.get_or("sliding_window", "vmin", "2")), "percentile_vmin");
    double clip_high_pct = finite_number(config.get_or("sliding_window", "percentile_vmax",
        config.get_or("sliding_window", "vmax", "98")), "percentile_vmax");
    if (clip_low_pct < 0 || clip_high_pct > 100 || clip_high_pct < 0 || clip_low_pct > 100)
        throw std::runtime_error("Clipping percentiles must lie in [0,100].");
    if (clip_low_pct > clip_high_pct) std::swap(clip_low_pct, clip_high_pct);
    const bool moment = statistic_name == "moment";
    const bool use_bell = evaluator.bell;
    const bool weighted_enabled = boolean(config.get_or("weighted_moment", "enable_weighted_moment",
        config.get_or("weighted_moment", "enabled", config.get_or("weighted_moment", "use_weighted_moment", "false")))) ||
        evaluator.weighted || source_weights.declustered;
    const bool center_cumulant = !moment && order >= 2 &&
        (center || source_weights.declustered || (weighted_enabled && options.zero_first_order_moment));
    const PointSpatialIndex index(data);
    std::vector<PointSample> locations;
    if (centers.regular) {
        locations.reserve(centers.count);
        for (size_t i=0;i<centers.count;++i) locations.push_back(centers.point<PointSample>(i));
    } else locations=data.nodes;
    std::vector<PointSample> output = locations;
    std::vector<PointSample> counts = locations;
    std::vector<PointSample> valid = locations;
    std::vector<std::string> reasons(locations.size());
    std::vector<double> finite_values;
    std::vector<WindowWeightSupport> weight_support(locations.size());

    for (size_t i = 0; i < locations.size(); ++i) {
        const PointSample& anchor = locations[i];
        const std::array<double, 3> position = {anchor.x, anchor.y, anchor.z};
        std::array<double, 3> lower, upper;
        for (size_t axis = 0; axis < 3; ++axis) {
            const double padding = 1e-12 * std::max(1.0, neighborhood.extent[axis]);
            lower[axis] = position[axis] - neighborhood.extent[axis] - padding;
            upper[axis] = position[axis] + neighborhood.extent[axis] + padding;
        }
        std::vector<double> samples;
        std::vector<double> weights;
        for (size_t candidate : index.query_box(lower, upper)) {
            // Identity-based exclusion preserves separate observations at duplicate coordinates.
            if (!centers.regular && !includes_center && candidate == i) continue;
            if (neighborhood.contains(data.nodes[candidate], anchor)) {
                samples.push_back(data.nodes[candidate].value);
                weights.push_back(source_weights.values[candidate]);
                weight_support[i].add(weights.back());
            }
        }
        counts[i].value = static_cast<double>(samples.size());
        counts[i].informed = true;
        output[i].value = std::numeric_limits<double>::quiet_NaN();
        output[i].informed = false;
        valid[i].value = 0.0;
        valid[i].informed = true;
        reasons[i] = samples.empty() ? "empty_neighborhood" : "minimum_samples";
        if (samples.size() < static_cast<size_t>(minimum)) continue;

        double result;
        if (moment) {
            result = weighted_moment_core::weighted_univariate_moment(samples.data(), weights.data(),
                static_cast<int>(samples.size()), order, true,
                weighted_enabled && options.zero_first_order_moment);
        }
        else {
            result = (use_bell ? bell : avl)(samples, order, center_cumulant,
                source_weights.declustered ? weights : std::vector<double>());
            if (order == 1 && weighted_enabled && options.zero_first_order_moment) result = 0.0;
        }
        if (invert) result = -result;
        if (!std::isfinite(result)) { reasons[i] = "nonfinite_output"; continue; }
        output[i].value = result;
        output[i].informed = true;
        valid[i].value = 1.0;
        reasons[i] = "valid";
        finite_values.push_back(result);
    }
    if (clip && !finite_values.empty()) {
        const double lo = percentile(finite_values, clip_low_pct), hi = percentile(finite_values, clip_high_pct);
        for (PointSample& sample : output)
            if (sample.informed) sample.value = std::clamp(sample.value, lo, hi);
    }
    std::filesystem::path outdir(config.get_or("runtime", "output_dir", "output/movewindow"));
    if (!outdir.is_absolute()) outdir = config_path.parent_path() / outdir;
    std::filesystem::create_directories(outdir);
    write_point_set(outdir / (target + ".txt"), data.title, target, output);
    write_point_set(outdir / (target + ".sample_count.txt"), data.title, target + "_sample_count", counts);
    write_point_set(outdir / (target + ".valid.txt"), data.title, target + "_valid", valid);
    source_weights.write(outdir, target, data);

    std::ofstream support(outdir / (target + ".window_support.csv"));
    std::ofstream meta(outdir / (target + ".meta.txt"));
    if (!support || !meta) throw std::runtime_error("Unable to write point-set window diagnostics.");
    support << (centers.regular ? "center_id,x,y,z,source_record_assigned," : "point_id,x,y,z,center_informed,")
            << "weight_sum,weight_squared_sum,effective_sample_size,"
               "sample_count,map_included,valid,valid_reason\n" << std::setprecision(17);
    for (size_t i = 0; i < locations.size(); ++i) {
        const auto& p = locations[i];
        const double n = counts[i].value;
        support << p.index << ',' << p.x << ',' << p.y << ',' << p.z << ',' << (p.informed ? 1 : 0)
            << ',' << static_cast<double>(weight_support[i].sum)
            << ',' << static_cast<double>(weight_support[i].squares)
            << ',' << weight_support[i].ess() << ',' << n << ",1," << valid[i].value << ',' << reasons[i] << '\n';
    }
    meta << std::setprecision(17)
         << "source_grid=" << source_path.string() << "\nsource_type=point_set\n"
         << "coordinate_dimensions=" << data.dimensions << "\ncoordinates=physical_unmodified\n"
         << "center_locations=" << (centers.regular ? "regular_grid" : "all_source_points")
         << "\nsource_points=" << data.nodes.size() << '\n'
         << "target_property=" << target << "\nstatistic=" << statistic << "\norder=" << order << '\n'
         << "evaluator=" << (moment ? "moment_average" : (use_bell ? "univariate_recurrence" : "AVL")) << '\n'
         << "weight_mode=" << template_weight_mode(config) << "\nnormalized_weights=true\n"
         << "point_weight_scope=source_dataset\nwindow_weight_normalization=selected_informed_samples\n"
         << source_weights.metadata
         << "weighted_moment_enabled=" << (weighted_enabled ? "true" : "false") << '\n'
         << "zero_first_order_moment=" << (options.zero_first_order_moment ? "true" : "false") << '\n'
         << "moment_centering=" << (moment && weighted_enabled && options.zero_first_order_moment ? "window_mean" : "none") << '\n'
         << "neighborhood_shape=" << (neighborhood.box ? "box" : "ellipsoid")
         << "\nneighborhood_units=source_coordinate_units\n"
         << "neighborhood_radii=" << neighborhood.radii[0] << ',' << neighborhood.radii[1] << ',' << neighborhood.radii[2] << '\n'
         << "neighborhood_angles=" << neighborhood.angles[0] << ',' << neighborhood.angles[1] << ',' << neighborhood.angles[2] << '\n'
         << "angle_convention=SGeMS_degrees_Rz(90-azimuth)_Ry(-dip)_Rx(rake)\n"
         << "includes_center=" << (includes_center ? "true" : "false") << "\nmin_samples=" << minimum << '\n'
         << "cumulant_centering=" << (center_cumulant ? (source_weights.declustered ? "weighted_window_mean" : "window_mean") : "none") << '\n'
         << "invalid_output=NaN\nwindow_support_schema=movewindow_point_support.v1\n"
         << "window_support_rows=" << output.size() << "\nvalid_output_rows=" << finite_values.size() << '\n'
         << "window_support_sidecar=" << target << ".window_support.csv\n"
         << "sample_count_sidecar=" << target << ".sample_count.txt\n"
         << "valid_sidecar=" << target << ".valid.txt\n";
    centers.metadata(meta);
    if (!support || !meta) throw std::runtime_error("Failed while writing point-set window diagnostics.");
    return 0;
}

} // namespace without_sgems
#endif
