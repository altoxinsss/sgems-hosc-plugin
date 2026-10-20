#include "config.h"
#include "grid_io.h"
#include "text.h"

#include "FunctionCall.h"
#include "Initialization.h"
#include "WeightedMomentCore.h"

#include <cstdlib>
#include <array>
#include <algorithm>
#include <cctype>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

namespace fs = std::filesystem;
using without_sgems::GridData;
using without_sgems::IniConfig;
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

bool has_alias(const IniConfig& config, const std::string& section, const std::vector<std::string>& keys) {
    for (size_t i = 0; i < keys.size(); ++i) {
        if (config.has(section, keys[i])) {
            return true;
        }
    }
    return false;
}

fs::path resolve_path(const fs::path& config_path, const std::string& value) {
    fs::path path(value);
    if (path.is_absolute()) {
        return path;
    }
    return config_path.parent_path() / path;
}

std::vector<int> parse_optional_int_list(const std::string& value) {
    if (without_sgems::trim(value).empty()) {
        return std::vector<int>();
    }
    return without_sgems::parse_int_list(value);
}

std::string join_ints(const std::vector<int>& values) {
    std::ostringstream stream;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i != 0) {
            stream << ",";
        }
        stream << values[i];
    }
    return stream.str();
}

void free_keys(Keys* key) {
    if (key == NULL) {
        return;
    }

    std::free(key->NOD);
    std::free(key->NLAG);
    std::free(key->LAG);
    std::free(key->LTOL);
    std::free(key->AZM);
    std::free(key->ATOL);
    std::free(key->BANDWH);
    std::free(key->DIP);
    std::free(key->DTOL);
    std::free(key->BANDWD);
    std::free(key);
}

GridData build_dense_grid_for_template(const GridData& source) {
    const size_t total = static_cast<size_t>(source.nx) * source.ny * source.nz;
    if (total == 0) {
        throw std::runtime_error("Source grid has an invalid dense extent.");
    }

    GridData dense = source;
    dense.nd = static_cast<int>(total);
    dense.nvar = 4;
    dense.nodes.clear();
    dense.nodes.reserve(total);

    std::vector<double> values(total, 0.0);
    for (size_t i = 0; i < source.nodes.size(); ++i) {
        const without_sgems::GridNode& node = source.nodes[i];
        const size_t index = static_cast<size_t>(node.z) * source.nx * source.ny +
            static_cast<size_t>(node.y) * source.nx +
            static_cast<size_t>(node.x);
        if (index >= values.size()) {
            throw std::runtime_error("Encountered an out-of-range node while densifying the template source grid.");
        }
        values[index] = node.value;
    }

    int index = 1;
    for (int z = 0; z < source.nz; ++z) {
        for (int y = 0; y < source.ny; ++y) {
            for (int x = 0; x < source.nx; ++x) {
                without_sgems::GridNode node;
                node.index = index++;
                node.x = x;
                node.y = y;
                node.z = z;
                node.value = values[static_cast<size_t>(z) * source.nx * source.ny +
                    static_cast<size_t>(y) * source.nx +
                    static_cast<size_t>(x)];
                dense.nodes.push_back(node);
            }
        }
    }

    return dense;
}

Keys* build_keys(
    const GridData& grid,
    const std::vector<int>& nod,
    const std::vector<int>& nlag,
    const std::vector<int>& lag,
    const std::vector<int>& ltol,
    const std::vector<int>& azm,
    const std::vector<int>& atol,
    const std::vector<int>& bandwh,
    const std::vector<int>& dip,
    const std::vector<int>& dtol,
    const std::vector<int>& bandwd) {
    const size_t dir_count = nod.size();
    if (dir_count == 0) {
        throw std::runtime_error("At least one direction must be configured for hosc_avl.");
    }

    const size_t expected = dir_count;
    const std::vector<const std::vector<int>*> lists = {
        &nlag, &lag, &ltol, &azm, &atol, &bandwh, &dip, &dtol, &bandwd
    };
    for (size_t i = 0; i < lists.size(); ++i) {
        if (lists[i]->size() != expected) {
            throw std::runtime_error("Directional parameter lists must all have the same length.");
        }
    }

    Keys* key = static_cast<Keys*>(std::calloc(1, sizeof(Keys)));
    if (key == NULL) {
        throw std::runtime_error("Failed to allocate Keys structure.");
    }

    key->DIM = (grid.nz <= 1) ? 2 : 3;
    key->REG_IRREG = 1;
    key->NX = grid.nx;
    key->NY = grid.ny;
    key->NZ = grid.nz;
    key->DIMX = (grid.nx > 0) ? grid.nx - 1 : 0;
    key->DIMY = (grid.ny > 0) ? grid.ny - 1 : 0;
    key->DIMZ = (grid.nz > 0) ? grid.nz - 1 : 0;
    key->ORDER = static_cast<int>(dir_count) + 1;
    key->NDIR = static_cast<int>(dir_count);

    key->NOD = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->NLAG = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->LAG = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->LTOL = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->AZM = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->ATOL = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->BANDWH = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->DIP = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->DTOL = static_cast<int*>(std::calloc(dir_count, sizeof(int)));
    key->BANDWD = static_cast<int*>(std::calloc(dir_count, sizeof(int)));

    if (key->NOD == NULL || key->NLAG == NULL || key->LAG == NULL || key->LTOL == NULL || key->AZM == NULL ||
        key->ATOL == NULL || key->BANDWH == NULL || key->DIP == NULL || key->DTOL == NULL || key->BANDWD == NULL) {
        free_keys(key);
        throw std::runtime_error("Failed to allocate directional arrays.");
    }

    for (size_t i = 0; i < dir_count; ++i) {
        key->NOD[i] = nod[i];
        key->NLAG[i] = nlag[i];
        key->LAG[i] = lag[i];
        key->LTOL[i] = ltol[i];
        key->AZM[i] = azm[i];
        key->ATOL[i] = atol[i];
        key->BANDWH[i] = bandwh[i];
        key->DIP[i] = dip[i];
        key->DTOL[i] = dtol[i];
        key->BANDWD[i] = bandwd[i];
    }

    return key;
}

struct CrossSectionLayout {
    int out_nx = 1;
    int out_ny = 1;
    int out_nz = 1;
    int dir_x = 0;
    int dir_y = 0;
    int dir_z = 0;
    std::vector<int> displayed_dirs;
    std::vector<int> hidden_dirs;
    std::vector<int> fixed_positions;
    std::vector<int> fixed_positions_zero_based;
};

CrossSectionLayout apply_cross_section_layout(
    Keys* key,
    const std::vector<int>& requested_displayed_dirs,
    const std::vector<int>& requested_positions) {
    if (key == NULL) {
        throw std::runtime_error("Keys must not be null.");
    }

    CrossSectionLayout layout;
    layout.displayed_dirs = requested_displayed_dirs;
    if (layout.displayed_dirs.size() != 2 && layout.displayed_dirs.size() != 3) {
        throw std::runtime_error("displayed_dirs must contain either two directions for a 2D section or three directions for a 3D section.");
    }

    std::vector<bool> displayed_mask(static_cast<size_t>(key->NDIR), false);
    for (size_t i = 0; i < layout.displayed_dirs.size(); ++i) {
        const int dir = layout.displayed_dirs[i];
        if (dir < 1 || dir > key->NDIR) {
            throw std::runtime_error("Displayed direction selection is out of range.");
        }
        if (displayed_mask[static_cast<size_t>(dir - 1)]) {
            throw std::runtime_error("Displayed direction selections must be distinct.");
        }
        displayed_mask[static_cast<size_t>(dir - 1)] = true;
    }

    for (int dir = 1; dir <= key->NDIR; ++dir) {
        if (!displayed_mask[static_cast<size_t>(dir - 1)]) {
            layout.hidden_dirs.push_back(dir);
        }
    }

    layout.fixed_positions = requested_positions;
    if (layout.fixed_positions.empty() && !layout.hidden_dirs.empty()) {
        layout.fixed_positions.assign(layout.hidden_dirs.size(), 1);
    }
    if (layout.hidden_dirs.empty()) {
        if (layout.fixed_positions.size() > 1 ||
            (layout.fixed_positions.size() == 1 && layout.fixed_positions[0] != 1)) {
            throw std::runtime_error("No fixed cross-section coordinates are needed when all directions are displayed.");
        }
        layout.fixed_positions.clear();
    }
    if (layout.fixed_positions.size() != layout.hidden_dirs.size()) {
        throw std::runtime_error("fixed_positions must contain one coordinate for each hidden direction.");
    }

    for (size_t i = 0; i < layout.fixed_positions.size(); ++i) {
        const int position = layout.fixed_positions[i];
        const int original_dir = layout.hidden_dirs[i];
        if (position < 1 || position > key->NLAG[original_dir - 1]) {
            throw std::runtime_error("Fixed cross-section coordinates are out of range.");
        }
        layout.fixed_positions_zero_based.push_back(position - 1);
    }

    layout.dir_x = layout.displayed_dirs[0];
    layout.dir_y = layout.displayed_dirs[1];
    layout.dir_z = (layout.displayed_dirs.size() == 3) ? layout.displayed_dirs[2] : 0;
    layout.out_nx = key->NLAG[layout.dir_x - 1];
    layout.out_ny = key->NLAG[layout.dir_y - 1];
    layout.out_nz = (layout.dir_z == 0) ? 1 : key->NLAG[layout.dir_z - 1];

    const std::vector<int> old_nod(key->NOD, key->NOD + key->NDIR);
    const std::vector<int> old_nlag(key->NLAG, key->NLAG + key->NDIR);
    const std::vector<int> old_lag(key->LAG, key->LAG + key->NDIR);
    const std::vector<int> old_azm(key->AZM, key->AZM + key->NDIR);
    const std::vector<int> old_atol(key->ATOL, key->ATOL + key->NDIR);
    const std::vector<int> old_bandwh(key->BANDWH, key->BANDWH + key->NDIR);
    const std::vector<int> old_dip(key->DIP, key->DIP + key->NDIR);
    const std::vector<int> old_dtol(key->DTOL, key->DTOL + key->NDIR);
    const std::vector<int> old_bandwd(key->BANDWD, key->BANDWD + key->NDIR);

    std::vector<int> reordered_dirs = layout.hidden_dirs;
    reordered_dirs.insert(reordered_dirs.end(), layout.displayed_dirs.begin(), layout.displayed_dirs.end());
    for (int i = 0; i < key->NDIR; ++i) {
        const int old_index = reordered_dirs[static_cast<size_t>(i)] - 1;
        key->NOD[i] = old_nod[static_cast<size_t>(old_index)];
        key->NLAG[i] = old_nlag[static_cast<size_t>(old_index)];
        key->LAG[i] = old_lag[static_cast<size_t>(old_index)];
        key->AZM[i] = old_azm[static_cast<size_t>(old_index)];
        key->ATOL[i] = old_atol[static_cast<size_t>(old_index)];
        key->BANDWH[i] = old_bandwh[static_cast<size_t>(old_index)];
        key->DIP[i] = old_dip[static_cast<size_t>(old_index)];
        key->DTOL[i] = old_dtol[static_cast<size_t>(old_index)];
        key->BANDWD[i] = old_bandwd[static_cast<size_t>(old_index)];
    }

    return layout;
}

bool is_moment_target(const std::string& value) {
    return without_sgems::to_lower(value).find("moment") != std::string::npos;
}

bool is_bell_algorithm(const std::string& value) {
    return without_sgems::to_lower(value).find("bell") != std::string::npos;
}

bool is_weighted_algorithm(const std::string& value) {
    return without_sgems::to_lower(value).find("weighted") != std::string::npos;
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

wmc::WeightMode parse_weight_mode(const std::string& value) {
    const std::string normalized = without_sgems::to_lower(value);
    if (normalized == "coordinate_gaussian" || normalized == "spatial_gaussian" || normalized == "coordinate") {
        return wmc::WeightMode::CoordinateGaussian;
    }
    if (normalized == "value_gaussian" || normalized == "value") {
        return wmc::WeightMode::ValueGaussian;
    }
    if (normalized == "soft_high_porosity" || normalized == "high_porosity") {
        return wmc::WeightMode::SoftHighPorosity;
    }
    if (normalized == "adaptive_upper_tail" || normalized == "adaptive_tail") {
        return wmc::WeightMode::AdaptiveUpperTail;
    }
    return wmc::WeightMode::Uniform;
}

wmc::WeightedMomentOptions parse_weight_options(const IniConfig& config) {
    wmc::WeightedMomentOptions options;
    const std::vector<std::string> scheme_switches = {
        "enable_coordinate_gaussian",
        "enable_fixed_value_gaussian",
        "enable_soft_high_porosity",
        "enable_adaptive_upper_tail"
    };
    bool has_explicit_switch = false;
    int enabled_scheme_count = 0;
    for (size_t i = 0; i < scheme_switches.size(); ++i) {
        if (config.has("weighted_moment", scheme_switches[i])) {
            has_explicit_switch = true;
            if (parse_bool_option(config.get("weighted_moment", scheme_switches[i]), false)) {
                ++enabled_scheme_count;
            }
        }
    }
    if (enabled_scheme_count > 1) {
        throw std::runtime_error(
            "Only one base weighting switch may be enabled: coordinate Gaussian, fixed value Gaussian, "
            "soft high porosity, or adaptive upper tail.");
    }
    if (has_explicit_switch) {
        if (parse_bool_option(get_alias_or(config, "weighted_moment", {"enable_coordinate_gaussian"}, "false"), false)) {
            options.mode = wmc::WeightMode::CoordinateGaussian;
        }
        else if (parse_bool_option(get_alias_or(config, "weighted_moment", {"enable_fixed_value_gaussian"}, "false"), false)) {
            options.mode = wmc::WeightMode::ValueGaussian;
        }
        else if (parse_bool_option(get_alias_or(config, "weighted_moment", {"enable_soft_high_porosity"}, "false"), false)) {
            options.mode = wmc::WeightMode::SoftHighPorosity;
        }
        else if (parse_bool_option(get_alias_or(config, "weighted_moment", {"enable_adaptive_upper_tail"}, "false"), false)) {
            options.mode = wmc::WeightMode::AdaptiveUpperTail;
        }
        else {
            options.mode = wmc::WeightMode::Uniform;
        }
    }
    else {
        options.mode = parse_weight_mode(get_alias_or(config, "weighted_moment", {"weight_mode"}, "uniform"));
    }
    options.bandwidth = std::atof(get_alias_or(config, "weighted_moment", {"bandwidth"}, "1").c_str());
    options.reference_value = std::atof(get_alias_or(config, "weighted_moment", {"reference_value"}, "0").c_str());
    options.min_weight = std::atof(get_alias_or(config, "weighted_moment", {"min_weight"}, "0").c_str());
    options.normalize_weights = parse_bool_option(
        get_alias_or(config, "weighted_moment", {"normalize_weights"}, "true"), true);
    options.zero_first_order_moment = parse_bool_option(
        get_alias_or(config, "weighted_moment", {"zero_first_order_moment"}, "false"), false);
    options.high_porosity_threshold = std::atof(get_alias_or(
        config, "weighted_moment", {"high_porosity_threshold"}, "0.25").c_str());
    options.soft_transition_width = std::atof(get_alias_or(
        config, "weighted_moment", {"soft_transition_width"}, "0.02").c_str());
    options.soft_weight_strength = std::atof(get_alias_or(
        config, "weighted_moment", {"soft_weight_strength"}, "1").c_str());
    options.adaptive_global_threshold = std::atof(get_alias_or(
        config, "weighted_moment", {"adaptive_global_threshold"}, "0.25").c_str());
    options.adaptive_mean_offset = std::atof(get_alias_or(
        config, "weighted_moment", {"adaptive_mean_offset"}, "1").c_str());
    options.adaptive_transition_scale = std::atof(get_alias_or(
        config, "weighted_moment", {"adaptive_transition_scale"}, "0.25").c_str());
    options.adaptive_weight_strength = std::atof(get_alias_or(
        config, "weighted_moment", {"adaptive_weight_strength"}, "1").c_str());
    options.adaptive_std_floor = std::atof(get_alias_or(
        config, "weighted_moment", {"adaptive_std_floor"}, "0.000001").c_str());
    options.joint_geometric_mean = parse_bool_option(get_alias_or(
        config, "weighted_moment", {"enable_joint_geometric_mean"}, "false"), false);

    if (options.mode == wmc::WeightMode::ValueGaussian && options.bandwidth <= 0.0) {
        throw std::runtime_error("bandwidth must be positive for fixed value Gaussian weighting.");
    }
    if (options.mode == wmc::WeightMode::SoftHighPorosity && options.soft_transition_width <= 0.0) {
        throw std::runtime_error("soft_transition_width must be positive.");
    }
    if (options.soft_weight_strength < 0.0 || options.soft_weight_strength > 1.0 ||
        options.adaptive_weight_strength < 0.0 || options.adaptive_weight_strength > 1.0) {
        throw std::runtime_error("soft_weight_strength and adaptive_weight_strength must lie in [0,1].");
    }
    if (options.mode == wmc::WeightMode::AdaptiveUpperTail &&
        (options.adaptive_transition_scale < 0.0 || options.adaptive_std_floor <= 0.0)) {
        throw std::runtime_error(
            "adaptive_transition_scale must be non-negative and adaptive_std_floor must be positive.");
    }
    if (options.joint_geometric_mean &&
        (options.mode == wmc::WeightMode::Uniform || options.mode == wmc::WeightMode::CoordinateGaussian)) {
        throw std::runtime_error(
            "enable_joint_geometric_mean is valid only for value-based weighting schemes.");
    }
    return options;
}

// The legacy HOSC search selects one unused node from an angular tolerance
// region.  Formal directional experiments need a fixed geometry instead: all
// valid translations of the same signed integer-offset template are retained.
// This implementation intentionally lives in template_cli so Bell and AVL use
// the same production template samples.
struct ExhaustiveSamplingLayout {
    std::vector<int> displayed_dirs;
    std::vector<int> hidden_dirs;
    std::vector<int> fixed_positions_zero_based;
};

struct ExhaustiveLagCell {
    std::vector<int> positions;
    std::vector<std::array<int, 3>> offsets;
};

struct ExhaustiveEvaluation {
    double value = std::numeric_limits<double>::quiet_NaN();
    int support_count = 0;
};

std::string normalized_template_sampling(const IniConfig& config) {
    return without_sgems::to_lower(get_alias_or(
        config, "template_shape", {"template_sampling", "sampling_mode"},
        get_alias_or(config, "general", {"template_sampling", "sampling_mode"}, "legacy_search")));
}

bool is_exhaustive_integer_sampling(const std::string& value) {
    return value == "exhaustive_integer_offsets" || value == "exhaustive_integer_offset";
}

std::vector<std::string> split_direction_vectors(const std::string& text) {
    std::vector<std::string> result;
    std::string current;
    for (size_t index = 0; index < text.size(); ++index) {
        const char ch = text[index];
        if (ch == ';' || ch == '|') {
            if (!without_sgems::trim(current).empty()) {
                result.push_back(without_sgems::trim(current));
            }
            current.clear();
        }
        else {
            current.push_back(ch);
        }
    }
    if (!without_sgems::trim(current).empty()) {
        result.push_back(without_sgems::trim(current));
    }
    return result;
}

std::vector<std::array<int, 3>> parse_integer_direction_vectors(
    const IniConfig& config,
    int coordinate_dimensions,
    size_t expected_count) {
    const std::string raw = get_alias(
        config, "template_shape", {"integer_direction_vectors", "direction_vectors", "signed_integer_directions"});
    const std::vector<std::string> tokens = split_direction_vectors(raw);
    if (tokens.size() != expected_count) {
        throw std::runtime_error("integer_direction_vectors must provide one vector for each configured direction.");
    }

    std::vector<std::array<int, 3>> result;
    result.reserve(tokens.size());
    for (size_t index = 0; index < tokens.size(); ++index) {
        std::string token = tokens[index];
        for (size_t character = 0; character < token.size(); ++character) {
            if (token[character] == '(' || token[character] == ')' || token[character] == '[' || token[character] == ']') {
                token[character] = ' ';
            }
        }
        const std::vector<int> parts = without_sgems::parse_int_list(token);
        if (parts.size() != static_cast<size_t>(coordinate_dimensions)) {
            throw std::runtime_error("Each integer direction vector must contain exactly coordinate_dimensions entries.");
        }
        std::array<int, 3> direction = {{0, 0, 0}};
        bool nonzero = false;
        for (int axis = 0; axis < coordinate_dimensions; ++axis) {
            direction[static_cast<size_t>(axis)] = parts[static_cast<size_t>(axis)];
            nonzero = nonzero || parts[static_cast<size_t>(axis)] != 0;
        }
        if (!nonzero) {
            throw std::runtime_error("integer_direction_vectors must not contain a zero vector.");
        }
        result.push_back(direction);
    }
    return result;
}

ExhaustiveSamplingLayout build_exhaustive_layout(
    const IniConfig& config,
    size_t direction_count,
    const std::vector<int>& nlag) {
    ExhaustiveSamplingLayout layout;
    if (has_alias(config, "target_cross_sections", {"displayed_dirs"})) {
        layout.displayed_dirs = parse_optional_int_list(config.get("target_cross_sections", "displayed_dirs"));
    }
    else {
        const int first = std::atoi(get_alias_or(config, "target_cross_sections", {"direction_x", "fst_dir"}, "0").c_str());
        const int second = std::atoi(get_alias_or(config, "target_cross_sections", {"direction_y", "snd_dir"}, "0").c_str());
        const int third = std::atoi(get_alias_or(config, "target_cross_sections", {"direction_z", "trd_dir"}, "0").c_str());
        if (first != 0) layout.displayed_dirs.push_back(first);
        if (second != 0) layout.displayed_dirs.push_back(second);
        if (third != 0) layout.displayed_dirs.push_back(third);
    }
    if (layout.displayed_dirs.empty() || layout.displayed_dirs.size() > 3) {
        throw std::runtime_error("exhaustive integer sampling requires one to three displayed directions.");
    }
    std::vector<bool> displayed(direction_count, false);
    for (size_t index = 0; index < layout.displayed_dirs.size(); ++index) {
        const int direction = layout.displayed_dirs[index];
        if (direction < 1 || direction > static_cast<int>(direction_count) || displayed[static_cast<size_t>(direction - 1)]) {
            throw std::runtime_error("displayed_dirs must be distinct, one-based configured direction indices.");
        }
        displayed[static_cast<size_t>(direction - 1)] = true;
    }
    for (size_t index = 0; index < direction_count; ++index) {
        if (!displayed[index]) layout.hidden_dirs.push_back(static_cast<int>(index) + 1);
    }

    std::vector<int> fixed = parse_optional_int_list(get_alias_or(
        config, "target_cross_sections", {"fixed_positions", "position", "cross_pos"}, ""));
    if (fixed.empty() && !layout.hidden_dirs.empty()) {
        fixed.assign(layout.hidden_dirs.size(), 1);
    }
    if (fixed.size() != layout.hidden_dirs.size()) {
        throw std::runtime_error("fixed_positions must contain one one-based lag position per hidden direction.");
    }
    for (size_t index = 0; index < fixed.size(); ++index) {
        const int direction = layout.hidden_dirs[index] - 1;
        if (fixed[index] < 1 || fixed[index] > nlag[static_cast<size_t>(direction)]) {
            throw std::runtime_error("fixed_positions contains an out-of-range lag position.");
        }
        layout.fixed_positions_zero_based.push_back(fixed[index] - 1);
    }
    return layout;
}

void enumerate_exhaustive_lag_cells(
    const std::vector<int>& nlag,
    const std::vector<int>& lag_step,
    const std::vector<std::array<int, 3>>& directions,
    const ExhaustiveSamplingLayout& layout,
    size_t direction_index,
    std::vector<int>* current,
    std::vector<ExhaustiveLagCell>* output) {
    if (current == NULL || output == NULL) {
        throw std::runtime_error("Internal exhaustive lag enumeration error.");
    }
    if (direction_index == nlag.size()) {
        ExhaustiveLagCell cell;
        cell.positions = *current;
        cell.offsets.reserve(directions.size());
        for (size_t index = 0; index < directions.size(); ++index) {
            const int distance = cell.positions[index] * lag_step[index];
            cell.offsets.push_back({{
                directions[index][0] * distance,
                directions[index][1] * distance,
                directions[index][2] * distance,
            }});
        }
        output->push_back(cell);
        return;
    }
    const std::vector<int>::const_iterator hidden = std::find(
        layout.hidden_dirs.begin(), layout.hidden_dirs.end(), static_cast<int>(direction_index) + 1);
    if (hidden != layout.hidden_dirs.end()) {
        const size_t hidden_index = static_cast<size_t>(hidden - layout.hidden_dirs.begin());
        (*current)[direction_index] = layout.fixed_positions_zero_based[hidden_index];
        enumerate_exhaustive_lag_cells(nlag, lag_step, directions, layout, direction_index + 1, current, output);
        return;
    }
    for (int position = 0; position < nlag[direction_index]; ++position) {
        (*current)[direction_index] = position;
        enumerate_exhaustive_lag_cells(nlag, lag_step, directions, layout, direction_index + 1, current, output);
    }
}

double dense_grid_value(const GridData& grid, int x, int y, int z) {
    const size_t index = static_cast<size_t>(z) * static_cast<size_t>(grid.nx) * static_cast<size_t>(grid.ny) +
        static_cast<size_t>(y) * static_cast<size_t>(grid.nx) + static_cast<size_t>(x);
    if (index >= grid.nodes.size()) {
        throw std::runtime_error("Exhaustive integer sampler encountered an out-of-range dense grid coordinate.");
    }
    return grid.nodes[index].value;
}

void free_exhaustive_templates(Template** templates, int template_count) {
    if (templates == NULL || *templates == NULL) return;
    Template* values = *templates;
    for (int index = 0; index < template_count; ++index) {
        std::free(values[index].Value);
        std::free(values[index].Lag);
    }
    std::free(values[0].index);
    std::free(values);
    *templates = NULL;
    G_Template = NULL;
    G_TemplateNum = 0;
}

Template* build_exhaustive_templates(
    const GridData& grid,
    const std::vector<std::array<int, 3>>& offsets,
    const std::vector<int>& nod,
    int* template_count) {
    if (template_count == NULL || offsets.size() != nod.size() || offsets.empty()) {
        throw std::runtime_error("Invalid exhaustive template geometry.");
    }
    const int variable_count = static_cast<int>(offsets.size()) + 1;
    std::array<int, 3> minimum = {{0, 0, 0}};
    std::array<int, 3> maximum = {{0, 0, 0}};
    for (size_t direction = 0; direction < offsets.size(); ++direction) {
        for (int axis = 0; axis < 3; ++axis) {
            minimum[static_cast<size_t>(axis)] = std::min(minimum[static_cast<size_t>(axis)], offsets[direction][static_cast<size_t>(axis)]);
            maximum[static_cast<size_t>(axis)] = std::max(maximum[static_cast<size_t>(axis)], offsets[direction][static_cast<size_t>(axis)]);
        }
    }
    const std::array<int, 3> dimensions = {{grid.nx, grid.ny, grid.nz}};
    for (int axis = 0; axis < 3; ++axis) {
        if (maximum[static_cast<size_t>(axis)] - minimum[static_cast<size_t>(axis)] >= dimensions[static_cast<size_t>(axis)]) {
            *template_count = 0;
            return NULL;
        }
    }

    const int capacity = (grid.nx - maximum[0] + minimum[0]) *
        (grid.ny - maximum[1] + minimum[1]) *
        (grid.nz - maximum[2] + minimum[2]);
    if (capacity <= 0) {
        *template_count = 0;
        return NULL;
    }
    Template* templates = static_cast<Template*>(std::calloc(static_cast<size_t>(capacity), sizeof(Template)));
    if (templates == NULL) throw std::runtime_error("Failed to allocate exhaustive integer templates.");

    int count = 0;
    for (int z = -minimum[2]; z < grid.nz - maximum[2]; ++z) {
        for (int y = -minimum[1]; y < grid.ny - maximum[1]; ++y) {
            for (int x = -minimum[0]; x < grid.nx - maximum[0]; ++x) {
                std::vector<double> row(static_cast<size_t>(variable_count), 0.0);
                row[0] = dense_grid_value(grid, x, y, z);
                bool finite = std::isfinite(row[0]);
                for (size_t direction = 0; direction < offsets.size() && finite; ++direction) {
                    const std::array<int, 3>& offset = offsets[direction];
                    const double value = dense_grid_value(grid, x + offset[0], y + offset[1], z + offset[2]);
                    finite = std::isfinite(value);
                    row[direction + 1] = value;
                }
                if (!finite) continue;
                templates[count].Value = static_cast<double*>(std::calloc(static_cast<size_t>(variable_count) + 1, sizeof(double)));
                templates[count].Lag = static_cast<int*>(std::calloc(static_cast<size_t>(variable_count) + 1, sizeof(int)));
                if (templates[count].Value == NULL || templates[count].Lag == NULL) {
                    free_exhaustive_templates(&templates, count + 1);
                    throw std::runtime_error("Failed to allocate an exhaustive integer template row.");
                }
                for (int variable = 0; variable < variable_count; ++variable) {
                    templates[count].Value[variable + 1] = row[static_cast<size_t>(variable)];
                }
                ++count;
            }
        }
    }
    if (count == 0) {
        std::free(templates);
        *template_count = 0;
        return NULL;
    }
    for (int variable = 1; variable <= variable_count; ++variable) {
        double mean = 0.0;
        for (int row = 0; row < count; ++row) mean += templates[row].Value[variable];
        mean /= static_cast<double>(count);
        for (int row = 0; row < count; ++row) templates[row].Value[variable] -= mean;
    }
    templates[0].index = static_cast<int*>(std::calloc(static_cast<size_t>(variable_count) + 1, sizeof(int)));
    if (templates[0].index == NULL) {
        free_exhaustive_templates(&templates, count);
        throw std::runtime_error("Failed to allocate exhaustive template multiplicities.");
    }
    int joint_order = 1;
    templates[0].index[1] = 1;
    for (size_t direction = 0; direction < nod.size(); ++direction) {
        if (nod[direction] <= 0) {
            free_exhaustive_templates(&templates, count);
            throw std::runtime_error("NOD entries must be positive for exhaustive integer sampling.");
        }
        templates[0].index[direction + 2] = nod[direction];
        joint_order += nod[direction];
    }
    templates[0].index[0] = joint_order;
    *template_count = count;
    return templates;
}

double evaluate_exhaustive_cumulant(Template* templates, int template_count, int variable_count, bool bell) {
    if (templates == NULL || template_count <= 0 || templates[0].index == NULL) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    G_Template = templates;
    G_TemplateNum = template_count;
    ResetMomentWorkspaceForDLL();
    ClearAVLWorkspaceForDLL(variable_count);
    double value = myEmpty;
    if (bell) {
        value = BellCumulantFromMultiplicityIndex(templates[0].index, variable_count);
    }
    else {
        SetFactorSwitchesForRun(1, 1, 1, 1);
        if (PrecomputeMomentsForIndex(templates[0].index, variable_count) != myEmpty &&
            InitAVLWorkspaceForDLL(variable_count) == 0) {
            value = BetterRecCumulant3(templates[0].index, variable_count);
        }
    }
    ClearOneEvaluationWorkspaceForDLL(variable_count);
    return value;
}

void write_exhaustive_template_output(
    const fs::path& output_dir,
    const GridData& source_grid,
    const std::string& target_property,
    const std::string& source_grid_path,
    const std::vector<double>& origin,
    const std::vector<int>& nlag,
    const std::vector<int>& lag_step,
    const std::vector<std::array<int, 3>>& directions,
    const ExhaustiveSamplingLayout& layout,
    const std::vector<ExhaustiveLagCell>& cells,
    const std::vector<ExhaustiveEvaluation>& results,
    const std::vector<double>& source_spacing,
    const std::string& algorithm) {
    if (cells.size() != results.size()) throw std::runtime_error("Exhaustive template output has inconsistent row counts.");
    fs::create_directories(output_dir);
    std::vector<without_sgems::GridNode> nodes;
    nodes.reserve(cells.size());
    for (size_t row = 0; row < cells.size(); ++row) {
        without_sgems::GridNode node;
        node.index = static_cast<int>(row) + 1;
        node.x = cells[row].positions[static_cast<size_t>(layout.displayed_dirs[0] - 1)];
        node.y = layout.displayed_dirs.size() > 1
            ? cells[row].positions[static_cast<size_t>(layout.displayed_dirs[1] - 1)] : 0;
        node.z = layout.displayed_dirs.size() > 2
            ? cells[row].positions[static_cast<size_t>(layout.displayed_dirs[2] - 1)] : 0;
        node.value = results[row].value;
        nodes.push_back(node);
    }
    without_sgems::write_text_grid(
        output_dir / (target_property + ".txt"), source_grid,
        source_grid.title.empty() ? target_property : source_grid.title, target_property, nodes);

    std::ofstream meta((output_dir / (target_property + ".meta.txt")).c_str(), std::ios::out | std::ios::trunc);
    meta << "source_grid=" << source_grid_path << "\n";
    meta << "target_property=" << target_property << "\n";
    meta << "target_origin=" << origin[0] << "," << origin[1] << "," << origin[2] << "\n";
    meta << "template_sampling=exhaustive_integer_offsets\n";
    meta << "calculate_algorithm=" << algorithm << "\n";
    meta << "displayed_dirs=" << join_ints(layout.displayed_dirs) << "\n";
    meta << "output_lag_counts=" << join_ints(nlag) << "\n";
    meta << "lag_steps_cells=" << join_ints(lag_step) << "\n";
    meta << "source_spacing=" << std::setprecision(17) << source_spacing[0] << "," << source_spacing[1] << "," << source_spacing[2] << "\n";
    meta << "integer_direction_vectors=";
    for (size_t direction = 0; direction < directions.size(); ++direction) {
        if (direction != 0) meta << ";";
        meta << directions[direction][0] << "," << directions[direction][1] << "," << directions[direction][2];
    }
    meta << "\n";

    std::ofstream support((output_dir / (target_property + ".support.csv")).c_str(), std::ios::out | std::ios::trunc);
    support << "lag_cell,lag_positions,integer_offsets,physical_offsets,valid_translations,support_count,cumulant\n";
    support << std::setprecision(17);
    for (size_t row = 0; row < cells.size(); ++row) {
        std::ostringstream positions, offsets, physical;
        for (size_t direction = 0; direction < directions.size(); ++direction) {
            if (direction != 0) {
                positions << ";";
                offsets << ";";
                physical << ";";
            }
            positions << cells[row].positions[direction];
            const std::array<int, 3>& offset = cells[row].offsets[direction];
            offsets << "(" << offset[0] << "," << offset[1] << "," << offset[2] << ")";
            physical << "(" << offset[0] * source_spacing[0] << "," << offset[1] * source_spacing[1] << "," << offset[2] * source_spacing[2] << ")";
        }
        support << row << ",\"" << positions.str() << "\",\"" << offsets.str() << "\",\"" << physical.str() << "\","
                << results[row].support_count << "," << results[row].support_count << "," << results[row].value << "\n";
    }
}

int run_exhaustive_integer_sampling(
    const IniConfig& config,
    const GridData& grid,
    const fs::path& config_path,
    const fs::path& source_grid_path,
    const std::string& calculate_target,
    const std::string& calculate_algorithm,
    const std::string& target_property,
    const std::vector<double>& origin,
    const std::vector<int>& nod,
    const std::vector<int>& nlag,
    const std::vector<int>& lag_step) {
    if (is_weighted_algorithm(calculate_algorithm) || is_moment_target(calculate_target)) {
        throw std::runtime_error("exhaustive_integer_offsets supports unweighted cumulant targets only.");
    }
    if (nod.empty() || nod.size() != nlag.size() || nod.size() != lag_step.size()) {
        throw std::runtime_error("NOD, NLAG, and LAG must have equal nonzero lengths for exhaustive integer sampling.");
    }
    const int coordinate_dimensions = grid.nz <= 1 ? 2 : 3;
    const std::vector<std::array<int, 3>> directions = parse_integer_direction_vectors(config, coordinate_dimensions, nod.size());
    const ExhaustiveSamplingLayout layout = build_exhaustive_layout(config, nod.size(), nlag);
    std::vector<ExhaustiveLagCell> cells;
    std::vector<int> positions(nod.size(), 0);
    enumerate_exhaustive_lag_cells(nlag, lag_step, directions, layout, 0, &positions, &cells);
    if (cells.empty()) throw std::runtime_error("No lag cells were selected for exhaustive integer sampling.");
    const bool bell = is_bell_algorithm(calculate_algorithm);
    std::vector<ExhaustiveEvaluation> results;
    results.reserve(cells.size());
    const int variable_count = static_cast<int>(nod.size()) + 1;
    for (size_t row = 0; row < cells.size(); ++row) {
        int support_count = 0;
        Template* templates = build_exhaustive_templates(grid, cells[row].offsets, nod, &support_count);
        if (templates == NULL || support_count <= 0) {
            throw std::runtime_error("An exhaustive lag cell has no complete template translations.");
        }
        const double value = evaluate_exhaustive_cumulant(templates, support_count, variable_count, bell);
        free_exhaustive_templates(&templates, support_count);
        if (value == myEmpty || !std::isfinite(value)) {
            throw std::runtime_error("Bell/AVL failed while evaluating an exhaustive integer template.");
        }
        ExhaustiveEvaluation evaluation;
        evaluation.value = value;
        evaluation.support_count = support_count;
        results.push_back(evaluation);
    }
    std::vector<double> source_spacing = without_sgems::parse_double_list(get_alias_or(
        config, "template_shape", {"source_spacing", "physical_spacing"},
        get_alias_or(config, "general", {"source_spacing", "physical_spacing"}, "1,1,1")));
    if (source_spacing.size() != 3 || source_spacing[0] <= 0.0 || source_spacing[1] <= 0.0 || source_spacing[2] <= 0.0) {
        throw std::runtime_error("source_spacing must contain three positive physical spacings.");
    }
    const fs::path output_dir = resolve_path(config_path, get_alias_or(config, "runtime", {"output_dir"}, "output/hosc_avl"));
    write_exhaustive_template_output(
        output_dir, grid, target_property, source_grid_path.string(), origin, nlag, lag_step, directions, layout,
        cells, results, source_spacing, calculate_algorithm);
    std::cout << "Wrote exhaustive integer-offset template output to " << output_dir.string() << "\n";
    return 0;
}

bool lag_matches_cross_section(const Keys* key, const int* lag, const int* crossection, int crossection_count) {
    if (key == NULL || lag == NULL) {
        return false;
    }
    if (crossection == NULL || crossection_count <= 0) {
        return true;
    }
    for (int i = 0; i < crossection_count; ++i) {
        if (crossection[i] < 0 || crossection[i] >= key->NLAG[i]) {
            return false;
        }
        if (lag[i] != crossection[i] * key->LAG[i]) {
            return false;
        }
    }
    return true;
}

double** save_selected_result_grid(
    const std::vector<double>& values,
    int** lags,
    const std::vector<int>& selected_lag_positions,
    int direction_count) {
    double** target_grid = static_cast<double**>(std::calloc(selected_lag_positions.size(), sizeof(double*)));
    if (target_grid == NULL) {
        return NULL;
    }
    for (size_t row = 0; row < selected_lag_positions.size(); ++row) {
        target_grid[row] = static_cast<double*>(std::calloc(static_cast<size_t>(direction_count) + 1, sizeof(double)));
        if (target_grid[row] == NULL) {
            for (size_t i = 0; i < row; ++i) {
                std::free(target_grid[i]);
            }
            std::free(target_grid);
            return NULL;
        }
        const int lag_position = selected_lag_positions[row];
        for (int j = 0; j < direction_count; ++j) {
            target_grid[row][j] = lags[lag_position][j];
        }
        target_grid[row][direction_count] = values[row];
    }
    return target_grid;
}

double calculate_weighted_template_statistic(
    Template* templates,
    int template_count,
    int variable_count,
    const std::vector<int>& multiplicities,
    const wmc::WeightedMomentOptions& options,
    bool moment_target) {
    if (templates == NULL || template_count <= 0 || variable_count <= 0) {
        return 0.0;
    }

    std::vector<std::vector<double>> sample_values(static_cast<size_t>(template_count));
    std::vector<wmc::WeightedSample> samples(static_cast<size_t>(template_count));
    for (int i = 0; i < template_count; ++i) {
        sample_values[static_cast<size_t>(i)].resize(static_cast<size_t>(variable_count), 0.0);
        for (int j = 0; j < variable_count; ++j) {
            sample_values[static_cast<size_t>(i)][static_cast<size_t>(j)] = templates[i].Value[j + 1];
        }
        samples[static_cast<size_t>(i)].values = sample_values[static_cast<size_t>(i)].data();
        samples[static_cast<size_t>(i)].value_count = variable_count;
        samples[static_cast<size_t>(i)].x = static_cast<double>(i);
        samples[static_cast<size_t>(i)].y = 0.0;
        samples[static_cast<size_t>(i)].z = 0.0;
    }

    if (moment_target) {
        return wmc::weighted_moment(
            samples.data(), template_count, multiplicities.data(), variable_count, options);
    }
    return wmc::weighted_cumulant_from_multiplicities(
        samples.data(), template_count, multiplicities.data(), variable_count, options);
}

double** WeightedMoment4Dll(
    Keys* key,
    double**& grid,
    int* output_rows,
    int* crossection,
    int crossection_count,
    const wmc::WeightedMomentOptions& options,
    bool moment_target);

int build_weighted_moment_key(const int* index, int order, char* out, int out_size) {
    if (index == NULL || out == NULL || out_size < order + 1) {
        return -1;
    }
    bool skipped = false;
    int zero_num = 0;
    for (int i = order; i >= 1; --i) {
        if (index[i] == 0 && skipped == false) {
            zero_num++;
            continue;
        }
        if (index[i] != 0 && skipped == false) {
            skipped = true;
        }
        out[i - 1] = static_cast<char>('0' + index[i]);
    }
    out[order - zero_num] = '\0';
    return 0;
}

double calculate_weighted_moment_for_index(
    Template* templates,
    int template_count,
    int variable_count,
    const int* index,
    const wmc::WeightedMomentOptions& options) {
    if (templates == NULL || template_count <= 0 || variable_count <= 0 || index == NULL) {
        return 0.0;
    }

    std::vector<int> multiplicities(static_cast<size_t>(variable_count), 0);
    for (int j = 0; j < variable_count; ++j) {
        multiplicities[static_cast<size_t>(j)] = index[j + 1];
    }

    std::vector<std::vector<double>> sample_values(static_cast<size_t>(template_count));
    std::vector<wmc::WeightedSample> samples(static_cast<size_t>(template_count));
    for (int i = 0; i < template_count; ++i) {
        sample_values[static_cast<size_t>(i)].resize(static_cast<size_t>(variable_count), 0.0);
        for (int j = 0; j < variable_count; ++j) {
            sample_values[static_cast<size_t>(i)][static_cast<size_t>(j)] = templates[i].Value[j + 1];
        }
        samples[static_cast<size_t>(i)].values = sample_values[static_cast<size_t>(i)].data();
        samples[static_cast<size_t>(i)].value_count = variable_count;
        samples[static_cast<size_t>(i)].x = static_cast<double>(i);
        samples[static_cast<size_t>(i)].y = 0.0;
        samples[static_cast<size_t>(i)].z = 0.0;
    }

    return wmc::weighted_moment(
        samples.data(), template_count, multiplicities.data(), variable_count, options);
}

int precompute_weighted_moments_for_avl(
    Template* templates,
    int template_count,
    const int* index,
    int order,
    const wmc::WeightedMomentOptions& options) {
    if (templates == NULL || template_count <= 0 || index == NULL || order <= 0) {
        return -1;
    }
    if (G_MomentFull == true) {
        return 0;
    }

    int* limits = static_cast<int*>(std::calloc(static_cast<size_t>(order) + 1, sizeof(int)));
    if (limits == NULL) {
        return -1;
    }

    int moment_count = 1;
    for (int j = 1; j <= order; ++j) {
        limits[j] = index[j];
        moment_count *= (limits[j] + 1);
    }

    G_Moment = static_cast<double*>(std::calloc(static_cast<size_t>(moment_count), sizeof(double)));
    G_MomentIndex = static_cast<char**>(std::calloc(static_cast<size_t>(moment_count), sizeof(char*)));
    if (G_Moment == NULL || G_MomentIndex == NULL) {
        std::free(limits);
        return -1;
    }
    G_SeriesNum = moment_count;
    G_CurrMomentPos = 0;

    G_Moment[0] = calculate_weighted_moment_for_index(templates, template_count, order, index, options);
    G_MomentIndex[0] = static_cast<char*>(std::calloc(static_cast<size_t>(order) + 1, sizeof(char)));
    if (G_MomentIndex[0] == NULL || build_weighted_moment_key(index, order, G_MomentIndex[0], order + 1) != 0) {
        std::free(limits);
        return -1;
    }
    G_CurrMomentPos = 1;

    if (moment_count > 1) {
        int radix = 1;
        int* base = static_cast<int*>(std::calloc(static_cast<size_t>(order) + 1, sizeof(int)));
        if (base == NULL) {
            std::free(limits);
            return -1;
        }
        for (int j = 1; j <= order; ++j) {
            base[j] = limits[j] + 1;
            radix *= base[j];
        }

        for (int code = 1; code < radix; ++code) {
            int* local_exp = static_cast<int*>(std::calloc(static_cast<size_t>(order) + 1, sizeof(int)));
            if (local_exp == NULL) {
                std::free(base);
                std::free(limits);
                return -1;
            }
            local_exp[0] = 0;
            int rem = code;
            for (int j = order; j >= 1; --j) {
                local_exp[j] = rem % base[j];
                rem /= base[j];
            }

            G_Moment[code] = calculate_weighted_moment_for_index(templates, template_count, order, local_exp, options);
            G_MomentIndex[code] = static_cast<char*>(std::calloc(static_cast<size_t>(order) + 1, sizeof(char)));
            if (G_MomentIndex[code] == NULL ||
                build_weighted_moment_key(local_exp, order, G_MomentIndex[code], order + 1) != 0) {
                std::free(local_exp);
                std::free(base);
                std::free(limits);
                return -1;
            }
            std::free(local_exp);
        }
        G_CurrMomentPos = radix;
        std::free(base);
    }

    std::free(limits);
    G_MomentFull = true;
    G_MomentHashReady = false;
    if (G_MomentHash != NULL) {
        std::free(G_MomentHash);
        G_MomentHash = NULL;
    }
    G_MomentHashSize = 0;
    return 0;
}

double** WeightedRecCalculation4DLL2(
    Keys* key,
    double**& grid,
    int* output_rows,
    int* crossection,
    int crossection_count,
    const wmc::WeightedMomentOptions& options,
    bool moment_target) {
    if (moment_target) {
        return WeightedMoment4Dll(key, grid, output_rows, crossection, crossection_count, options, true);
    }
    if (key == NULL || grid == NULL || output_rows == NULL) {
        return NULL;
    }
    if (InitialStruct4(key, grid) != 0) {
        return NULL;
    }

    int** lags = NULL;
    int lag_num = 0;
    if (InitialLagMat(&lags, &lag_num) == -1) {
        Destruction();
        return NULL;
    }

    std::vector<int> selected_lag_positions;
    for (int i = 0; i < lag_num; ++i) {
        if (lag_matches_cross_section(key, lags[i], crossection, crossection_count)) {
            selected_lag_positions.push_back(i);
        }
    }
    if (*output_rows > 0 && static_cast<int>(selected_lag_positions.size()) != *output_rows) {
        Free2DVec_int(&lags, lag_num);
        Destruction();
        return NULL;
    }

    std::vector<double> values(selected_lag_positions.size(), 0.0);
    for (int row = static_cast<int>(selected_lag_positions.size()) - 1; row >= 0; --row) {
        Template* templates = NULL;
        int template_count = 0;
        const int lag_position = selected_lag_positions[static_cast<size_t>(row)];
        const int search_result = (key->DIM == 2)
            ? TemplateConstruct_2D(lags[lag_position], key->ORDER, &templates, &template_count, 0, pt_Grids->ND)
            : TemplateConstruct_3D(lags[lag_position], key->ORDER, &templates, &template_count, 0, pt_Grids->ND);
        if (search_result != 0) {
            continue;
        }

        ResetMomentWorkspaceForDLL();
        ClearAVLWorkspaceForDLL(key->ORDER);
        if (precompute_weighted_moments_for_avl(templates, template_count, templates->index, key->ORDER, options) != 0 ||
            InitAVLWorkspaceForDLL(key->ORDER) != 0) {
            FreeTemplate(&templates);
            Free2DVec_int(&lags, lag_num);
            Destruction();
            return NULL;
        }

        SetLowOrderCumulantMomentCacheMode(1);
        values[static_cast<size_t>(row)] = BetterRecCumulant3(templates->index, key->ORDER);
        SetLowOrderCumulantMomentCacheMode(0);
        FreeTemplate(&templates);
        ClearOneEvaluationWorkspaceForDLL(key->ORDER);

        if (values[static_cast<size_t>(row)] == myEmpty) {
            Free2DVec_int(&lags, lag_num);
            Destruction();
            return NULL;
        }
    }

    *output_rows = static_cast<int>(selected_lag_positions.size());
    double** result = save_selected_result_grid(values, lags, selected_lag_positions, key->NDIR);
    Free2DVec_int(&lags, lag_num);
    Destruction();
    return result;
}

double** WeightedMoment4Dll(
    Keys* key,
    double**& grid,
    int* output_rows,
    int* crossection,
    int crossection_count,
    const wmc::WeightedMomentOptions& options,
    bool moment_target) {
    if (key == NULL || grid == NULL || output_rows == NULL) {
        return NULL;
    }
    if (InitialStruct4(key, grid) != 0) {
        return NULL;
    }

    int** lags = NULL;
    int lag_num = 0;
    if (InitialLagMat(&lags, &lag_num) == -1) {
        Destruction();
        return NULL;
    }

    std::vector<int> selected_lag_positions;
    for (int i = 0; i < lag_num; ++i) {
        if (lag_matches_cross_section(key, lags[i], crossection, crossection_count)) {
            selected_lag_positions.push_back(i);
        }
    }
    if (*output_rows > 0 && static_cast<int>(selected_lag_positions.size()) != *output_rows) {
        Free2DVec_int(&lags, lag_num);
        Destruction();
        return NULL;
    }

    std::vector<int> multiplicities(static_cast<size_t>(key->ORDER), 0);
    multiplicities[0] = 1;
    for (int i = 1; i < key->ORDER; ++i) {
        multiplicities[static_cast<size_t>(i)] = key->NOD[i - 1];
    }

    std::vector<double> values(selected_lag_positions.size(), 0.0);
    for (int row = static_cast<int>(selected_lag_positions.size()) - 1; row >= 0; --row) {
        Template* templates = NULL;
        int template_count = 0;
        const int lag_position = selected_lag_positions[static_cast<size_t>(row)];
        const int search_result = (key->DIM == 2)
            ? TemplateConstruct_2D(lags[lag_position], key->ORDER, &templates, &template_count, 0, pt_Grids->ND)
            : TemplateConstruct_3D(lags[lag_position], key->ORDER, &templates, &template_count, 0, pt_Grids->ND);
        if (search_result != 0) {
            continue;
        }
        values[static_cast<size_t>(row)] = calculate_weighted_template_statistic(
            templates, template_count, key->ORDER, multiplicities, options, moment_target);
        FreeTemplate(&templates);
    }

    *output_rows = static_cast<int>(selected_lag_positions.size());
    double** result = save_selected_result_grid(values, lags, selected_lag_positions, key->NDIR);
    Free2DVec_int(&lags, lag_num);
    Destruction();
    return result;
}

void write_template_output(
    const fs::path& output_dir,
    const GridData& source_grid_style,
    const std::string& target_property,
    double** result,
    int rows,
    int direction_count,
    const CrossSectionLayout& layout,
    const std::vector<double>& origin,
    const fs::path& source_grid) {
    fs::create_directories(output_dir);

    const int hidden_dim_count = static_cast<int>(layout.hidden_dirs.size());
    const int x_column = hidden_dim_count;
    const int y_column = hidden_dim_count + 1;
    const int z_column = (layout.displayed_dirs.size() == 3) ? hidden_dim_count + 2 : -1;

    std::vector<without_sgems::GridNode> nodes;
    nodes.reserve(static_cast<size_t>(rows));
    for (int row = 0; row < rows; ++row) {
        without_sgems::GridNode node;
        node.index = row + 1;
        node.x = static_cast<int>(result[row][x_column]);
        node.y = static_cast<int>(result[row][y_column]);
        node.z = (z_column >= 0) ? static_cast<int>(result[row][z_column]) : 0;
        node.value = result[row][direction_count];
        nodes.push_back(node);
    }
    const fs::path txt_path = output_dir / (target_property + ".txt");
    const std::string output_title = source_grid_style.title.empty() ? target_property : source_grid_style.title;
    without_sgems::write_text_grid(txt_path, source_grid_style, output_title, target_property, nodes);

    const fs::path meta_path = output_dir / (target_property + ".meta.txt");
    std::ofstream meta(meta_path.c_str(), std::ios::out | std::ios::trunc);
    meta << "source_grid=" << source_grid.string() << "\n";
    meta << "target_property=" << target_property << "\n";
    meta << "target_origin=" << origin[0] << "," << origin[1] << "," << origin[2] << "\n";
    meta << "display_direction_x=" << layout.dir_x << "\n";
    meta << "display_direction_y=" << layout.dir_y << "\n";
    meta << "display_direction_z=" << layout.dir_z << "\n";
    meta << "displayed_dirs=" << join_ints(layout.displayed_dirs) << "\n";
    meta << "hidden_dirs=" << join_ints(layout.hidden_dirs) << "\n";
    meta << "fixed_positions=" << join_ints(layout.fixed_positions) << "\n";
    meta << "output_shape=" << layout.out_nx << "," << layout.out_ny << "," << layout.out_nz << "\n";
}

}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: hosc_avl_cli <config.ini>\n";
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

        std::vector<double> origin;
        if (config.has("general", "target_origin")) {
            origin = without_sgems::parse_double_list(config.get("general", "target_origin"));
        }
        else {
            origin.push_back(std::atof(get_alias_or(config, "general", {"ox", "xmn_"}, "0").c_str()));
            origin.push_back(std::atof(get_alias_or(config, "general", {"oy", "ymn_"}, "0").c_str()));
            origin.push_back(std::atof(get_alias_or(config, "general", {"oz", "zmn_"}, "0").c_str()));
        }
        if (origin.size() != 3) {
            throw std::runtime_error("target_origin must contain exactly three values.");
        }

        const GridData source_grid = without_sgems::load_hosc_grid(
            source_grid_path, source_property, coordinate_dimensions);
        const GridData grid = build_dense_grid_for_template(source_grid);
        const std::vector<int> nod = without_sgems::parse_int_list(get_alias(
            config, "template_shape", {"number_of_occurrences_in_each_direction", "nod"}));
        const std::vector<int> nlag = without_sgems::parse_int_list(get_alias(
            config, "template_shape", {"number_of_lags_in_each_direction", "nlag"}));
        const std::vector<int> lag = without_sgems::parse_int_list(get_alias(
            config, "template_shape", {"lag_separation_in_each_direction", "lag"}));
        const std::vector<int> ltol = without_sgems::parse_int_list(get_alias(
            config, "template_shape", {"lag_tolerance_in_each_direction", "ltol"}));
        const std::vector<int> azm = without_sgems::parse_int_list(get_alias(
            config, "template_shape", {"azimuth", "azm"}));
        const std::vector<int> atol = without_sgems::parse_int_list(get_alias(
            config, "template_shape", {"tolerance_of_azimuth", "atol"}));
        const std::vector<int> bandwh = without_sgems::parse_int_list(get_alias(
            config, "template_shape", {"bandwidth_of_azimuth", "bandwh"}));
        const std::vector<int> dip = without_sgems::parse_int_list(get_alias(
            config, "template_shape", {"dip_angle", "dip"}));
        const std::vector<int> dtol = without_sgems::parse_int_list(get_alias(
            config, "template_shape", {"tolerance_of_dip_angle", "dtol"}));
        const std::vector<int> bandwd = without_sgems::parse_int_list(get_alias(
            config, "template_shape", {"dip_bandwidth", "bandwd"}));

        const std::string sampling = normalized_template_sampling(config);
        if (is_exhaustive_integer_sampling(sampling)) {
            return run_exhaustive_integer_sampling(
                config, grid, config_path, source_grid_path, calculate_target, calculate_algorithm,
                target_property, origin, nod, nlag, lag);
        }
        if (sampling != "legacy_search" && sampling != "legacy" && sampling != "search") {
            throw std::runtime_error(
                "template_sampling must be legacy_search or exhaustive_integer_offsets.");
        }

        Keys* key = build_keys(grid, nod, nlag, lag, ltol, azm, atol, bandwh, dip, dtol, bandwd);

        try {
            std::vector<int> displayed_dirs;
            if (has_alias(config, "target_cross_sections", {"displayed_dirs"})) {
                displayed_dirs = parse_optional_int_list(config.get("target_cross_sections", "displayed_dirs"));
            }
            else {
                const int dir_x = std::atoi(get_alias(config, "target_cross_sections", {"direction_x", "fst_dir"}).c_str());
                const int dir_y = std::atoi(get_alias(config, "target_cross_sections", {"direction_y", "snd_dir"}).c_str());
                const int dir_z = std::atoi(get_alias_or(config, "target_cross_sections", {"direction_z", "trd_dir"}, "0").c_str());
                if (dir_x != 0) {
                    displayed_dirs.push_back(dir_x);
                }
                if (dir_y != 0) {
                    displayed_dirs.push_back(dir_y);
                }
                if (dir_z != 0) {
                    displayed_dirs.push_back(dir_z);
                }
            }
            const std::string positions_text = get_alias_or(
                config, "target_cross_sections", {"fixed_positions", "position", "cross_pos"}, "");
            const std::vector<int> positions = parse_optional_int_list(positions_text);
            const CrossSectionLayout layout = apply_cross_section_layout(key, displayed_dirs, positions);

            const fs::path output_dir = resolve_path(
                config_path, get_alias_or(config, "runtime", {"output_dir"}, "output/hosc_avl"));
            int output_rows = layout.out_nx * layout.out_ny * layout.out_nz;
            double** vendor_grid = without_sgems::create_vendor_grid_array(grid);
            double** result = NULL;
            std::vector<int> avl_cross_values = layout.fixed_positions_zero_based;

            const bool weighted_calculation = is_weighted_algorithm(calculate_algorithm);
            const wmc::WeightedMomentOptions weighted_options = parse_weight_options(config);
            if (weighted_calculation) {
                if (is_moment_target(calculate_target)) {
                    result = WeightedMoment4Dll(
                        key, vendor_grid, &output_rows, avl_cross_values.empty() ? NULL : avl_cross_values.data(),
                        static_cast<int>(avl_cross_values.size()), weighted_options, true);
                }
                else if (is_bell_algorithm(calculate_algorithm)) {
                    result = WeightedMoment4Dll(
                        key, vendor_grid, &output_rows, avl_cross_values.empty() ? NULL : avl_cross_values.data(),
                        static_cast<int>(avl_cross_values.size()), weighted_options, false);
                }
                else {
                    result = WeightedRecCalculation4DLL2(
                        key, vendor_grid, &output_rows, avl_cross_values.empty() ? NULL : avl_cross_values.data(),
                        static_cast<int>(avl_cross_values.size()), weighted_options, false);
                }
            }
            else if (is_moment_target(calculate_target)) {
                result = Moment4Dll(
                    key,
                    vendor_grid,
                    &output_rows,
                    avl_cross_values.empty() ? NULL : avl_cross_values.data(),
                    static_cast<int>(avl_cross_values.size()));
            }
            else if (is_bell_algorithm(calculate_algorithm)) {
                result = RecCalculation4DLL3(
                    key,
                    vendor_grid,
                    &output_rows,
                    avl_cross_values.empty() ? NULL : avl_cross_values.data(),
                    static_cast<int>(avl_cross_values.size()));
            }
            else {
                result = RecCalculation4DLL2(
                    key,
                    vendor_grid,
                    &output_rows,
                    avl_cross_values.empty() ? NULL : avl_cross_values.data(),
                    static_cast<int>(avl_cross_values.size()));
            }

            if (result == NULL) {
                throw std::runtime_error("The core hosc_avl calculation returned no output.");
            }

            write_template_output(output_dir, grid, target_property, result, output_rows, key->NDIR, layout, origin, source_grid_path);
            if (weighted_calculation) {
                const fs::path metadata_path = output_dir / (target_property + ".meta.txt");
                std::ofstream metadata(metadata_path.c_str(), std::ios::out | std::ios::app);
                metadata << "weighted_moment_enabled=true\n";
                metadata << wmc::weight_options_manifest(weighted_options);
            }
            without_sgems::free_vendor_matrix(result, output_rows);
            free_keys(key);

            std::cout << "Wrote hosc_avl output to " << output_dir.string() << "\n";
            return 0;
        }
        catch (...) {
            free_keys(key);
            throw;
        }
    }
    catch (const std::exception& error) {
        std::cerr << "hosc_avl_cli failed: " << error.what() << "\n";
        return 1;
    }
}
