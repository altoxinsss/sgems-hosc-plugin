#include "config.h"
#include "dense_grid.h"
#include "grid_io.h"
#include "text.h"
#include "point_window.h"
#include "weight_config.h"

#include "FunctionCall.h"
#include "Initialization.h"
#include "Tools.h"
#include "WeightedMomentCore.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using without_sgems::GridData;
using without_sgems::IniConfig;
using without_sgems::DenseGrid;
namespace wmc = weighted_moment_core;

namespace {

std::string get_alias(const IniConfig& config, const std::string& section, const std::vector<std::string>& keys) {
    for (size_t i = 0; i < keys.size(); ++i) {
        if (config.has(section, keys[i])) {
            return config.get(section, keys[i]);
        }
    }
    throw std::runtime_error("Missing configuration value in section [" + section + "]");
}

std::string get_alias_or(
    const IniConfig& config,
    const std::string& section,
    const std::vector<std::string>& keys,
    const std::string& fallback) {
    for (size_t i = 0; i < keys.size(); ++i) {
        if (config.has(section, keys[i])) {
            return config.get(section, keys[i]);
        }
    }
    return fallback;
}

fs::path resolve_path(const fs::path& config_path, const std::string& value) {
    fs::path path(value);
    if (path.is_absolute()) {
        return path;
    }
    return config_path.parent_path() / path;
}

bool is_moment_target(const std::string& value) {
    return without_sgems::to_lower(value).find("moment") != std::string::npos;
}

bool is_bell_algorithm(const std::string& value) {
    return without_sgems::to_lower(value).find("bell") != std::string::npos;
}

bool is_weighted_algorithm(const std::string& value) {
    const std::string normalized = without_sgems::to_lower(value);
    return normalized.find("weighted") != std::string::npos;
}

bool parse_bool_option(const std::string& value, bool fallback) {
    const std::string normalized = without_sgems::to_lower(value);
    if (normalized == "1" || normalized == "true" || normalized == "yes" || normalized == "on") {
        return true;
    }
    if (normalized == "0" || normalized == "false" || normalized == "no" || normalized == "off") {
        return false;
    }
    return fallback;
}

double clamp_double_value(double value, double lower, double upper) {
    if (value < lower) {
        return lower;
    }
    if (value > upper) {
        return upper;
    }
    return value;
}

double percentile_value(std::vector<double> values, double percentile) {
    if (values.empty()) {
        return 0.0;
    }
    percentile = std::clamp(percentile, 0.0, 100.0);
    std::sort(values.begin(), values.end());
    const double position = percentile / 100.0 * static_cast<double>(values.size() - 1);
    const size_t lower_index = static_cast<size_t>(std::floor(position));
    const size_t upper_index = static_cast<size_t>(std::ceil(position));
    if (lower_index == upper_index) {
        return values[lower_index];
    }
    const double fraction = position - static_cast<double>(lower_index);
    return values[lower_index] * (1.0 - fraction) + values[upper_index] * fraction;
}

bool is_weighted_enabled(const IniConfig& config) {
    if (config.has("weighted_moment", "enable_weighted_moment")) {
        return parse_bool_option(config.get("weighted_moment", "enable_weighted_moment"), false);
    }
    if (config.has("weighted_moment", "enabled")) {
        return parse_bool_option(config.get("weighted_moment", "enabled"), false);
    }
    if (config.has("weighted_moment", "use_weighted_moment")) {
        return parse_bool_option(config.get("weighted_moment", "use_weighted_moment"), false);
    }
    return false;
}

wmc::WeightedMomentOptions parse_weight_options(const IniConfig& config) {
    return without_sgems::parse_checked_weight_options(config, true);
}

double bell_binomial(int n, int k) {
    if (k < 0 || k > n) {
        return 0.0;
    }
    if (k == 0 || k == n) {
        return 1.0;
    }
    if (k > n - k) {
        k = n - k;
    }

    long double result = 1.0;
    for (int i = 1; i <= k; ++i) {
        result *= static_cast<long double>(n - k + i);
        result /= static_cast<long double>(i);
    }
    return static_cast<double>(result);
}

std::vector<double> centered_cumulant_samples(
    const std::vector<double>& samples,
    int order,
    bool enable_centering,
    const std::vector<double>& weights) {
    if (!enable_centering || order < 2 || samples.empty()) {
        return samples;
    }

    long double sum = 0.0, total = 0.0;
    for (size_t i = 0; i < samples.size(); ++i) {
        const long double weight = weights.empty() ? 1.0 : weights.at(i);
        sum += weight * samples[i];
        total += weight;
    }
    const double mean = static_cast<double>(sum / total);
    std::vector<double> centered;
    centered.reserve(samples.size());
    for (const double sample : samples) {
        centered.push_back(sample - mean);
    }
    return centered;
}

double calculate_window_bell_cumulant(
    const std::vector<double>& samples,
    int order,
    bool enable_centering,
    const std::vector<double>& weights) {
    if (order <= 0 || samples.empty()) {
        return 0.0;
    }

    const std::vector<double> stable_samples = centered_cumulant_samples(
        samples, order, enable_centering, weights);
    if (!weights.empty()) {
        return wmc::weighted_univariate_cumulant(stable_samples.data(), weights.data(),
            static_cast<int>(samples.size()), order, true, false);
    }

    std::vector<double> moments(static_cast<size_t>(order) + 1, 0.0);
    std::vector<double> cumulants(static_cast<size_t>(order) + 1, 0.0);

    for (size_t i = 0; i < stable_samples.size(); ++i) {
        long double power = 1.0;
        for (int n = 1; n <= order; ++n) {
            power *= static_cast<long double>(stable_samples[i]);
            moments[static_cast<size_t>(n)] += static_cast<double>(power);
        }
    }

    for (int n = 1; n <= order; ++n) {
        moments[static_cast<size_t>(n)] /= static_cast<double>(stable_samples.size());
    }

    cumulants[1] = moments[1];
    for (int n = 2; n <= order; ++n) {
        long double value = moments[static_cast<size_t>(n)];
        for (int j = 1; j < n; ++j) {
            value -= static_cast<long double>(bell_binomial(n - 1, j - 1)) *
                static_cast<long double>(cumulants[static_cast<size_t>(j)]) *
                static_cast<long double>(moments[static_cast<size_t>(n - j)]);
        }
        cumulants[static_cast<size_t>(n)] = static_cast<double>(value);
    }

    return cumulants[static_cast<size_t>(order)];
}

void clear_avl_roots(int order) {
    if (RootAVL2 == NULL) {
        return;
    }
    for (int i = 1; i <= order; ++i) {
        if (RootAVL2[i] != NULL) {
            FreeAVLTree(&RootAVL2[i]);
            RootAVL2[i] = NULL;
        }
    }
    std::free(RootAVL2);
    RootAVL2 = NULL;
}

void reset_moment_workspace() {
    if (G_Moment != NULL) {
        std::free(G_Moment);
        G_Moment = NULL;
    }
    if (G_MomentIndex != NULL) {
        if (G_SeriesNum > 0) {
            Free2DVec_char(&G_MomentIndex, G_SeriesNum);
        }
        else {
            std::free(G_MomentIndex);
            G_MomentIndex = NULL;
        }
    }
    if (G_MomentHash != NULL) {
        std::free(G_MomentHash);
        G_MomentHash = NULL;
    }
    G_MomentHashSize = 0;
    G_MomentHashReady = false;
    G_MomentFull = false;
    G_SeriesNum = 0;
    G_CurrMomentPos = 0;
}

double calculate_window_avl_cumulant(
    const std::vector<double>& samples,
    int order,
    bool enable_centering,
    const std::vector<double>& weights) {
    if (samples.empty()) {
        return 0.0;
    }

    const std::vector<double> stable_samples = centered_cumulant_samples(
        samples, order, enable_centering, weights);
    reset_moment_workspace();
    clear_avl_roots(1);
    struct Cleanup {
        ~Cleanup() {
            clear_avl_roots(1);
            reset_moment_workspace();
            SetLowOrderCumulantMomentCacheMode(0);
            ResetFactorSwitchesForRun();
            G_Template = NULL;
            G_TemplateNum = 0;
        }
    } cleanup;
    std::vector<std::array<double, 2>> values(samples.size());
    std::vector<Template> views(samples.size());
    int index[] = {order, order};
    for (size_t i = 0; i < samples.size(); ++i) {
        values[i] = {{0.0, stable_samples[i]}};
        views[i].Value = values[i].data();
        views[i].index = index;
    }
    G_Template = views.data();
    G_TemplateNum = static_cast<int>(views.size());
    if (!weights.empty()) {
        // Supply every moment needed by AVL, including its low-order base
        // cases. Slot zero is the root moment, as required by the evaluator.
        SetFactorSwitchesForRun(1, 1, 1, 1);
        SetLowOrderCumulantMomentCacheMode(1);
        G_SeriesNum = order + 1;
        G_Moment = static_cast<double*>(std::calloc(G_SeriesNum, sizeof(double)));
        G_MomentIndex = static_cast<char**>(std::calloc(G_SeriesNum, sizeof(char*)));
        if (!G_Moment || !G_MomentIndex) throw std::bad_alloc();
        for (int slot = 0; slot <= order; ++slot) {
            const int power = slot == 0 ? order : slot - 1;
            G_Moment[slot] = wmc::weighted_univariate_moment(stable_samples.data(),
                weights.data(), G_TemplateNum, power, true, false);
            if (!std::isfinite(G_Moment[slot]))
                throw std::runtime_error("Weighted window moment overflow.");
            G_MomentIndex[slot] = static_cast<char*>(std::calloc(2, sizeof(char)));
            if (!G_MomentIndex[slot]) throw std::bad_alloc();
            if (power) G_MomentIndex[slot][0] = static_cast<char>('0' + power);
        }
        G_CurrMomentPos = G_SeriesNum;
        G_MomentFull = true;
    }
    if (InitRoot2(1) == -1) {
        throw std::runtime_error("InitRoot2 failed for moving-window AVL mode.");
    }

    const double result = BetterRecCumulant3(G_Template->index, 1);
    RefreshTree(1);
    return result;
}

std::vector<double> parse_origin(const IniConfig& config) {
    if (config.has("general", "target_origin")) {
        return without_sgems::parse_double_list(config.get("general", "target_origin"));
    }

    std::vector<double> origin;
    origin.push_back(std::atof(get_alias_or(config, "general", {"ox", "xmn_"}, "0").c_str()));
    origin.push_back(std::atof(get_alias_or(config, "general", {"oy", "ymn_"}, "0").c_str()));
    origin.push_back(std::atof(get_alias_or(config, "general", {"oz", "zmn_"}, "0").c_str()));
    return origin;
}

struct WindowSupportRecord {
    int x = 0;
    int y = 0;
    int z = 0;
    double weight_sum = 0.0;
    double weight_squared_sum = 0.0;
    double effective_sample_size = 0.0;
    int sample_count = 0;
    bool map_included = false;
    bool valid = false;
    std::string valid_reason;
};

const char* window_support_schema_version() {
    return "movewindow_window_support.v1";
}

void write_window_support_sidecar(
    const fs::path& output_dir,
    const std::string& target_property,
    const std::vector<WindowSupportRecord>& records) {
    const fs::path sidecar_path = output_dir / (target_property + ".window_support.csv");
    std::ofstream sidecar(sidecar_path.c_str(), std::ios::out | std::ios::trunc);
    if (!sidecar) {
        throw std::runtime_error("Unable to create moving-window support sidecar: " + sidecar_path.string());
    }
    sidecar << "x,y,z,weight_sum,weight_squared_sum,effective_sample_size,"
            << "sample_count,map_included,valid,valid_reason\n";
    sidecar << std::setprecision(17);
    for (const WindowSupportRecord& record : records) {
        sidecar << record.x << ',' << record.y << ',' << record.z << ','
                << record.weight_sum << ',' << record.weight_squared_sum << ','
                << record.effective_sample_size << ',' << record.sample_count << ','
                << (record.map_included ? 1 : 0) << ',' << (record.valid ? 1 : 0) << ','
                << record.valid_reason << '\n';
    }
}

void write_window_output(
    const fs::path& output_dir,
    const GridData& source_grid_style,
    const std::string& target_property,
    const std::vector<double>& values,
    int out_nx,
    int out_ny,
    int out_nz,
    int step_x,
    int step_y,
    int step_z,
    int coordinate_offset_x,
    int coordinate_offset_y,
    int coordinate_offset_z,
    const std::vector<unsigned char>& keep_mask,
    const std::vector<double>& origin,
    const fs::path& source_grid) {
    fs::create_directories(output_dir);

    std::vector<without_sgems::GridNode> nodes;
    nodes.reserve(values.size());
    int block = 0;
    for (int z = 0; z < out_nz; ++z) {
        for (int y = 0; y < out_ny; ++y) {
            for (int x = 0; x < out_nx; ++x) {
                without_sgems::GridNode node;
                node.index = block + 1;
                node.x = source_grid_style.raw_min_x + x * step_x + coordinate_offset_x;
                node.y = source_grid_style.raw_min_y + y * step_y + coordinate_offset_y;
                node.z = source_grid_style.raw_min_z + z * step_z + coordinate_offset_z;
                node.value = values[static_cast<size_t>(block)];
                if (keep_mask.empty() || keep_mask[static_cast<size_t>(block)] != 0) {
                    nodes.push_back(node);
                }
                block++;
            }
        }
    }
    const fs::path txt_path = output_dir / (target_property + ".txt");
    const std::string output_title = source_grid_style.title.empty() ? target_property : source_grid_style.title;
    without_sgems::write_text_grid(txt_path, source_grid_style, output_title, target_property, nodes);

    const fs::path meta_path = output_dir / (target_property + ".meta.txt");
    std::ofstream meta(meta_path.c_str(), std::ios::out | std::ios::trunc);
    meta << "source_grid=" << source_grid.string() << "\n";
    meta << "target_property=" << target_property << "\n";
    meta << "target_origin=" << origin[0] << "," << origin[1] << "," << origin[2] << "\n";
    meta << "output_shape=" << out_nx << "," << out_ny << "," << out_nz << "\n";
    meta << "step=" << step_x << "," << step_y << "," << step_z << "\n";
    meta << "coordinate_offset=" << coordinate_offset_x << "," << coordinate_offset_y << "," << coordinate_offset_z << "\n";
    meta << "sparse_output_nodes=" << nodes.size() << "\n";
}

}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: movewindow_cli <config.ini>\n";
        return 1;
    }

    const fs::path config_path = fs::absolute(argv[1]);

    try {
        const IniConfig config = IniConfig::load(config_path);
        const std::string calculate_target = get_alias(config, "general", {"calculate_target", "target_button"});
        const std::string calculate_algorithm = get_alias_or(
            config, "general", {"calculate_algorithm", "algorithm_button"}, "AVL Based algorithm");
        const fs::path source_grid_path = resolve_path(
            config_path, get_alias(config, "general", {"source_grid", "source_grid_file"}));
        const std::string source_property = get_alias(config, "general", {"source_property"});
        const int coordinate_dimensions = std::atoi(get_alias_or(
            config, "general", {"coordinate_dimensions"}, "0").c_str());
        const std::string target_property = get_alias(config, "general", {"target_property", "props"});

        // Reject removed schemes before either grid loader can reinterpret data.
        const wmc::WeightedMomentOptions weighted_options = parse_weight_options(config);
        if (without_sgems::uses_point_set(config)) {
            return without_sgems::run_point_window(
                config, config_path, source_grid_path, source_property, coordinate_dimensions,
                target_property, calculate_target, calculate_algorithm,
                calculate_window_bell_cumulant, calculate_window_avl_cumulant);
        }

        const GridData sparse_grid = without_sgems::load_hosc_grid(
            source_grid_path, source_property, coordinate_dimensions, without_sgems::external_weight_property(config));
        const DenseGrid grid = without_sgems::build_dense_grid(sparse_grid);
        const without_sgems::DenseGridView grid_view(grid);
        const std::vector<double> source_origin = config.has("general", "source_origin")
            ? without_sgems::parse_double_list(config.get("general", "source_origin"))
            : std::vector<double>{
                  static_cast<double>(sparse_grid.raw_min_x),
                  static_cast<double>(sparse_grid.raw_min_y),
                  static_cast<double>(sparse_grid.raw_min_z)};
        const std::vector<double> source_spacing = config.has("general", "source_spacing")
            ? without_sgems::parse_double_list(config.get("general", "source_spacing"))
            : std::vector<double>{1.0, 1.0, 1.0};
        if (source_origin.size() != 3 || source_spacing.size() != 3 ||
            source_spacing[0] <= 0.0 || source_spacing[1] <= 0.0 || source_spacing[2] <= 0.0) {
            throw std::runtime_error("source_origin and source_spacing must each contain three values; spacing must be positive.");
        }

        const std::vector<int> window_size = config.has("sliding_window", "window_size")
            ? without_sgems::parse_int_list(config.get("sliding_window", "window_size"))
            : std::vector<int>{
                  std::atoi(get_alias(config, "sliding_window", {"windowl"}).c_str()),
                  std::atoi(get_alias(config, "sliding_window", {"windoww"}).c_str()),
                  std::atoi(get_alias_or(config, "sliding_window", {"windowh"}, "1").c_str())
              };
        const std::vector<int> move_step = config.has("sliding_window", "move_step")
            ? without_sgems::parse_int_list(config.get("sliding_window", "move_step"))
            : std::vector<int>{
                  std::atoi(get_alias(config, "sliding_window", {"stepx"}).c_str()),
                  std::atoi(get_alias(config, "sliding_window", {"stepy"}).c_str()),
                  std::atoi(get_alias_or(config, "sliding_window", {"stepz"}, "0").c_str())
              };
        if (window_size.size() != 3 || move_step.size() != 3) {
            throw std::runtime_error("window_size and move_step must each contain three values.");
        }

        const int order = std::atoi(get_alias(config, "sliding_window", {"moment_cumulant_order", "moment_order"}).c_str());
        if (order <= 0) {
            throw std::runtime_error("Moment/Cumulant order must be positive.");
        }
        const bool center_cumulant_input_for_order_ge_2 = parse_bool_option(
            get_alias_or(
                config,
                "sliding_window",
                {"center_cumulant_input_for_order_ge_2"},
                "false"),
            false);
        const bool invert_window_output = parse_bool_option(
            get_alias_or(config, "sliding_window", {"invert_window_output"}, "false"), false);
        const bool enable_percentile_clipping = parse_bool_option(
            get_alias_or(config, "sliding_window", {"enable_percentile_clipping", "percentile_clipping"}, "false"), false);
        double percentile_vmin = std::atof(
            get_alias_or(config, "sliding_window", {"percentile_vmin", "vmin"}, "2").c_str());
        double percentile_vmax = std::atof(
            get_alias_or(config, "sliding_window", {"percentile_vmax", "vmax"}, "98").c_str());
        const bool require_full_window = parse_bool_option(
            get_alias_or(config, "sliding_window", {"require_full_window"}, "true"), true);
        const bool require_complete_window = parse_bool_option(
            get_alias_or(config, "sliding_window", {"require_complete_window"}, "true"), true);
        const bool skip_missing_samples = parse_bool_option(
            get_alias_or(config, "sliding_window", {"skip_missing_samples"}, "false"), false);
        const bool require_anchor_present = parse_bool_option(
            get_alias_or(config, "sliding_window", {"require_anchor_present"}, "true"), true);
        const bool output_at_window_center = parse_bool_option(
            get_alias_or(config, "sliding_window", {"output_at_window_center"}, "true"), true);
        const bool write_window_support = parse_bool_option(
            get_alias_or(config, "runtime", {"write_window_support_sidecar"}, "true"), true);
        const int minimum_samples = std::max(
            1, std::atoi(get_alias_or(config, "sliding_window", {"minimum_samples"}, "1").c_str()));
        if (percentile_vmin > percentile_vmax) {
            std::swap(percentile_vmin, percentile_vmax);
        }

        int window_x = window_size[0];
        int window_y = window_size[1];
        int window_z = window_size[2];
        int step_x = move_step[0];
        int step_y = move_step[1];
        int step_z = move_step[2];

        if (step_x <= 0 || step_y <= 0) {
            throw std::runtime_error("Move steps in X and Y must be positive.");
        }

        const bool is_2d = (grid.nz <= 1);
        if (is_2d) {
            if (window_z != 0 && window_z != 1) {
                throw std::runtime_error("Window Z must be 0 or 1 for 2D grids.");
            }
            if (step_z != 0) {
                throw std::runtime_error("Move Step Z must be 0 for 2D grids.");
            }
            window_z = 1;
        }
        else if (step_z <= 0) {
            throw std::runtime_error("Move Step Z must be positive for 3D grids.");
        }

        if (require_full_window &&
            (window_x > grid.nx || window_y > grid.ny || (!is_2d && window_z > grid.nz))) {
            throw std::runtime_error("A full moving window does not fit inside the source grid.");
        }

        const int out_nx = require_full_window
            ? 1 + (grid.nx - window_x) / step_x
            : 1 + (grid.nx - 1) / step_x;
        const int out_ny = require_full_window
            ? 1 + (grid.ny - window_y) / step_y
            : 1 + (grid.ny - 1) / step_y;
        const int out_nz = is_2d
            ? 1
            : (require_full_window
                ? 1 + (grid.nz - window_z) / step_z
                : 1 + (grid.nz - 1) / step_z);

        const int coordinate_offset_x = output_at_window_center ? (window_x - 1) / 2 : 0;
        const int coordinate_offset_y = output_at_window_center ? (window_y - 1) / 2 : 0;
        const int coordinate_offset_z = (output_at_window_center && !is_2d) ? (window_z - 1) / 2 : 0;
        std::vector<double> output(static_cast<size_t>(out_nx) * out_ny * out_nz, 0.0);
        std::vector<unsigned char> keep_output(output.size(), 0);
        std::vector<WindowSupportRecord> window_support;
        window_support.reserve(output.size());
        std::vector<double> weight_sums;
        std::vector<double> weight_square_sums;
        std::vector<double> effective_sample_sizes;
        const bool declustered = without_sgems::uses_source_weights(config);
        without_sgems::PointSetData physical;
        physical.dimensions = is_2d ? 2 : 3;
        if (declustered) {
            const double missing = config.has("general", "missing_value")
                ? without_sgems::point_window_detail::finite_number(config.get("general", "missing_value"), "missing_value")
                : std::numeric_limits<double>::quiet_NaN();
            for (int z = 0; z < grid.nz; ++z) for (int y = 0; y < grid.ny; ++y) for (int x = 0; x < grid.nx; ++x) {
                without_sgems::PointSample p;
                p.index = grid.source_ids[physical.nodes.size()];
                p.external_weight = grid.external_weights[physical.nodes.size()];
                p.x = source_origin[0] + x * source_spacing[0];
                p.y = source_origin[1] + y * source_spacing[1];
                p.z = source_origin[2] + z * source_spacing[2];
                p.value = grid_view.value(x, y, z);
                p.informed = grid_view.present(x, y, z) && std::isfinite(p.value) && p.value != missing;
                physical.nodes.push_back(p);
            }
        }
        const without_sgems::WindowWeights source_weights(config, physical);
        const bool use_weighted_moment = declustered || is_weighted_enabled(config) || is_weighted_algorithm(calculate_algorithm);
        const bool cumulant_centering_applied =
            !is_moment_target(calculate_target) &&
            (use_weighted_moment
                ? (weighted_options.zero_first_order_moment ||
                    ((declustered || center_cumulant_input_for_order_ge_2) && order >= 2))
                : (center_cumulant_input_for_order_ge_2 && order >= 2));
        int block = 0;

        for (int anchor_z = 0; anchor_z < out_nz; ++anchor_z) {
            const int start_z = anchor_z * (is_2d ? 0 : step_z);
            for (int anchor_y = 0; anchor_y < out_ny; ++anchor_y) {
                const int start_y = anchor_y * step_y;
                for (int anchor_x = 0; anchor_x < out_nx; ++anchor_x) {
                    const int start_x = anchor_x * step_x;
                    const bool anchor_present = grid_view.present(start_x, start_y, start_z);
                    keep_output[static_cast<size_t>(block)] =
                        (!require_anchor_present || anchor_present) ? 1 : 0;
                    WindowSupportRecord support_record;
                    support_record.x = sparse_grid.raw_min_x + start_x + coordinate_offset_x;
                    support_record.y = sparse_grid.raw_min_y + start_y + coordinate_offset_y;
                    support_record.z = sparse_grid.raw_min_z + start_z + coordinate_offset_z;
                    support_record.map_included = keep_output[static_cast<size_t>(block)] != 0;
                    support_record.valid_reason = anchor_present
                        ? "valid"
                        : (support_record.map_included ? "anchor_missing_allowed" : "anchor_missing");
                    const int length = std::min(window_x, grid.nx - start_x);
                    const int width = std::min(window_y, grid.ny - start_y);
                    const int height = is_2d ? 1 : std::min(window_z, grid.nz - start_z);

                    std::vector<double> samples;
                    std::vector<double> weights;
                    bool complete_data_window = true;
                    samples.reserve(static_cast<size_t>(length) * width * height);
                    weights.reserve(static_cast<size_t>(length) * width * height);
                    for (int z = 0; z < height; ++z) {
                        for (int y = 0; y < width; ++y) {
                            for (int x = 0; x < length; ++x) {
                                const int gx = start_x + x;
                                const int gy = start_y + y;
                                const int gz = start_z + z;
                                const size_t source_id = (static_cast<size_t>(gz) * grid.ny + gy) * grid.nx + gx;
                                if (!grid_view.present(gx, gy, gz) ||
                                    (declustered && !physical.nodes[source_id].informed)) {
                                    complete_data_window = false;
                                    if (declustered || require_complete_window || skip_missing_samples) {
                                        continue;
                                    }
                                }
                                const double value = grid_view.value(gx, gy, gz);
                                samples.push_back(value);
                                if (use_weighted_moment) weights.push_back(declustered ? source_weights.values[source_id] : 1.0);
                            }
                        }
                    }

                    support_record.sample_count = static_cast<int>(samples.size());

                    if ((require_complete_window && !complete_data_window) ||
                        samples.size() < static_cast<size_t>(minimum_samples)) {
                        keep_output[static_cast<size_t>(block)] = 0;
                        output[static_cast<size_t>(block)] = 0.0;
                        support_record.map_included = false;
                        support_record.valid = false;
                        support_record.valid_reason = !complete_data_window
                            ? "incomplete_window"
                            : "minimum_samples";
                        window_support.push_back(support_record);
                        block++;
                        continue;
                    }

                    if (use_weighted_moment) {
                        long double weight_sum = 0.0;
                        long double weight_square_sum = 0.0;
                        for (const double weight : weights) {
                            weight_sum += weight;
                            weight_square_sum += weight * weight;
                        }
                        support_record.weight_sum = static_cast<double>(weight_sum);
                        support_record.weight_squared_sum = static_cast<double>(weight_square_sum);
                        if (weight_sum > 0.0 && weight_square_sum > 0.0) {
                            support_record.effective_sample_size = static_cast<double>(
                                weight_sum * weight_sum / weight_square_sum);
                            weight_sums.push_back(static_cast<double>(weight_sum));
                            weight_square_sums.push_back(static_cast<double>(weight_square_sum));
                            effective_sample_sizes.push_back(support_record.effective_sample_size);
                        }
                        if (is_moment_target(calculate_target)) {
                            output[static_cast<size_t>(block)] = wmc::weighted_univariate_moment(
                                samples.data(), weights.data(), static_cast<int>(samples.size()), order,
                                weighted_options.normalize_weights, weighted_options.zero_first_order_moment);
                        }
                        else {
                            const auto evaluator = without_sgems::point_window_detail::parse_evaluator(calculate_algorithm);
                            output[static_cast<size_t>(block)] = (evaluator.bell
                                ? calculate_window_bell_cumulant : calculate_window_avl_cumulant)(
                                    samples, order, cumulant_centering_applied, weights);
                            if (order == 1 && weighted_options.zero_first_order_moment)
                                output[static_cast<size_t>(block)] = 0.0;
                        }
                    }
                    else if (is_moment_target(calculate_target)) {
                        support_record.weight_sum = static_cast<double>(samples.size());
                        support_record.weight_squared_sum = static_cast<double>(samples.size());
                        support_record.effective_sample_size = static_cast<double>(samples.size());
                        long double sum = 0.0;
                        for (size_t i = 0; i < samples.size(); ++i) {
                            sum += std::pow(samples[i], static_cast<double>(order));
                        }
                        output[static_cast<size_t>(block)] = static_cast<double>(sum / samples.size());
                    }
                    else if (is_bell_algorithm(calculate_algorithm)) {
                        support_record.weight_sum = static_cast<double>(samples.size());
                        support_record.weight_squared_sum = static_cast<double>(samples.size());
                        support_record.effective_sample_size = static_cast<double>(samples.size());
                        output[static_cast<size_t>(block)] = calculate_window_bell_cumulant(
                            samples, order, cumulant_centering_applied, {});
                    }
                    else {
                        support_record.weight_sum = static_cast<double>(samples.size());
                        support_record.weight_squared_sum = static_cast<double>(samples.size());
                        support_record.effective_sample_size = static_cast<double>(samples.size());
                        output[static_cast<size_t>(block)] = calculate_window_avl_cumulant(
                            samples, order, cumulant_centering_applied, {});
                    }

                    if (invert_window_output) {
                        output[static_cast<size_t>(block)] = -output[static_cast<size_t>(block)];
                    }

                    if (!complete_data_window && support_record.valid_reason == "valid") {
                        support_record.valid_reason = "incomplete_window_allowed";
                    }
                    if (support_record.weight_sum <= 0.0 || support_record.weight_squared_sum <= 0.0) {
                        support_record.valid_reason = "nonpositive_weight_sum";
                    }
                    else if (!std::isfinite(output[static_cast<size_t>(block)])) {
                        support_record.valid_reason = "nonfinite_output";
                    }
                    support_record.valid = support_record.map_included && support_record.valid_reason == "valid";
                    window_support.push_back(support_record);

                    block++;
                }
            }
        }

        if (enable_percentile_clipping && !output.empty()) {
            std::vector<double> clipping_values;
            clipping_values.reserve(output.size());
            for (size_t i = 0; i < output.size(); ++i) {
                if (keep_output.empty() || keep_output[i] != 0) {
                    clipping_values.push_back(output[i]);
                }
            }
            if (clipping_values.empty()) {
                clipping_values = output;
            }
            const double lower_value = percentile_value(clipping_values, percentile_vmin);
            const double upper_value = percentile_value(clipping_values, percentile_vmax);
            const double clip_lower = std::min(lower_value, upper_value);
            const double clip_upper = std::max(lower_value, upper_value);
            for (size_t i = 0; i < output.size(); ++i) {
                output[i] = clamp_double_value(output[i], clip_lower, clip_upper);
            }
        }

        const fs::path output_dir = resolve_path(
            config_path, get_alias_or(config, "runtime", {"output_dir"}, "output/movewindow"));
        std::vector<double> output_origin = source_origin;
        if (output_at_window_center) {
            output_origin[0] += 0.5 * static_cast<double>(window_x - 1) * source_spacing[0];
            output_origin[1] += 0.5 * static_cast<double>(window_y - 1) * source_spacing[1];
            output_origin[2] += 0.5 * static_cast<double>((is_2d ? 1 : window_z) - 1) * source_spacing[2];
        }
        write_window_output(
            output_dir, sparse_grid, target_property, output, out_nx, out_ny, out_nz,
            step_x, step_y, step_z,
            coordinate_offset_x, coordinate_offset_y, coordinate_offset_z,
            keep_output, output_origin, source_grid_path);
        if (write_window_support) {
            write_window_support_sidecar(output_dir, target_property, window_support);
        }
        source_weights.write(output_dir, target_property, physical);
        {
            const fs::path metadata_path = output_dir / (target_property + ".meta.txt");
            std::ofstream metadata(metadata_path.c_str(), std::ios::out | std::ios::app);
            if (!metadata) {
                throw std::runtime_error("Unable to append moving-window metadata: " + metadata_path.string());
            }
            metadata << "weighted_moment_enabled=" << (use_weighted_moment ? "true" : "false") << "\n";
            metadata << wmc::weight_options_manifest(weighted_options);
            metadata << "window_weight_mode=" << without_sgems::template_weight_mode(config) << '\n'
                     << "point_weight_scope=source_dataset\nwindow_weight_normalization=selected_informed_samples\n"
                     << "evaluator=" << (is_moment_target(calculate_target) ? "moment_average" :
                        (without_sgems::point_window_detail::parse_evaluator(calculate_algorithm).bell ? "univariate_recurrence" : "AVL")) << '\n'
                     << source_weights.metadata;
            metadata << "center_cumulant_input_for_order_ge_2="
                     << (center_cumulant_input_for_order_ge_2 ? "true" : "false") << "\n";
            metadata << "cumulant_centering="
                     << (cumulant_centering_applied
                             ? (use_weighted_moment ? "weighted_window_mean" : "window_mean")
                             : "none")
                     << "\n";
            metadata << "cumulant_centering_version=order_ge_2.v1\n";
            metadata << "window_support_schema=" << window_support_schema_version() << "\n";
            metadata << "window_support_rows=" << window_support.size() << "\n";
            metadata << "window_support_sidecar="
                     << (write_window_support ? target_property + ".window_support.csv" : "disabled") << "\n";
            if (use_weighted_moment) {
            const auto write_diagnostic_summary = [&metadata](
                const char* name, const std::vector<double>& values) {
                if (values.empty()) {
                    return;
                }
                const auto range = std::minmax_element(values.begin(), values.end());
                const double sum = std::accumulate(values.begin(), values.end(), 0.0);
                metadata << name << "_min=" << *range.first << "\n";
                metadata << name << "_mean=" << sum / static_cast<double>(values.size()) << "\n";
                metadata << name << "_max=" << *range.second << "\n";
            };
            write_diagnostic_summary("weight_sum", weight_sums);
            write_diagnostic_summary("weight_square_sum", weight_square_sums);
            write_diagnostic_summary("effective_sample_size", effective_sample_sizes);
            }
        }
        std::cout << "Wrote movewindow output to " << output_dir.string() << "\n";
        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "movewindow_cli failed: " << error.what() << "\n";
        return 1;
    }
}
