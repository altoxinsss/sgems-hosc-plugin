#ifndef WITHOUT_SGEMS_WEIGHT_CONFIG_H
#define WITHOUT_SGEMS_WEIGHT_CONFIG_H

#include "config.h"
#include "text.h"
#include "WeightedMomentCore.h"

#include <stdexcept>
#include <string>

namespace without_sgems {

inline bool parse_weight_boolean(const std::string& value, const std::string& key) {
    const std::string normalized = to_lower(trim(value));
    if (normalized == "true" || normalized == "1" || normalized == "yes" || normalized == "on") {
        return true;
    }
    if (normalized == "false" || normalized == "0" || normalized == "no" || normalized == "off") {
        return false;
    }
    throw std::runtime_error("Invalid Boolean value for weighted_moment." + key + ": " + value);
}

// Always call this parser, even when weighted_moment.enabled is false, so
// obsolete mode names cannot silently execute as a different statistical run.
inline weighted_moment_core::WeightedMomentOptions parse_checked_weight_options(
    const IniConfig& config, bool allow_declustering = false) {
    const std::string section = "weighted_moment";
    const std::string mode = to_lower(trim(config.get_or(section, "weight_mode", "uniform")));
    const bool declustering = mode == "cell_declustering" || mode == "polygonal_declustering" || mode == "kriging_declustering";
    const bool external = mode == "external";
    if (mode != "uniform" && !(allow_declustering && (declustering || external))) {
        if (declustering) {
            throw std::runtime_error(mode + " is unavailable in this uniform-only calculation.");
        }
        if (mode == "external") {
            throw std::runtime_error(
                "External weights are unavailable in this uniform-only calculation.");
        }
        throw std::runtime_error("Unsupported weight_mode '" + mode +
            "'. Use uniform, cell_declustering, polygonal_declustering, kriging_declustering or external.");
    }

    const char* removed_switches[] = {
        "enable_coordinate_gaussian", "enable_fixed_value_gaussian",
        "enable_value_gaussian", "enable_soft_high_porosity",
        "enable_adaptive_upper_tail", "enable_joint_geometric_mean"
    };
    for (const char* key : removed_switches) {
        if (config.has(section, key) && parse_weight_boolean(config.get(section, key), key)) {
            throw std::runtime_error(std::string("Removed custom weighting switch: ") + key +
                ". Use uniform calculation or a supported declustering configuration.");
        }
    }

    weighted_moment_core::WeightedMomentOptions options;
    if (external && (trim(config.get_or(section, "weight_property", "")).empty() ||
                     trim(config.get_or("general", "source_property", "")).empty()))
        throw std::runtime_error("External weights require weighted_moment.weight_property and general.source_property.");
    if (declustering || external) {
        options.mode = weighted_moment_core::WeightMode::External;
        const char* enabled_keys[] = {"enabled", "enable_weighted_moment", "use_weighted_moment"};
        for (const char* key : enabled_keys) {
            if (config.has(section, key) && !parse_weight_boolean(config.get(section, key), key)) {
                throw std::runtime_error(mode + " conflicts with a disabled weighted_moment switch.");
            }
        }
        if (config.has("general", "enable_weighted_moment") &&
            !parse_weight_boolean(config.get("general", "enable_weighted_moment"), "enable_weighted_moment")) {
            throw std::runtime_error(mode + " conflicts with general.enable_weighted_moment=false.");
        }
    }
    if (config.has(section, "normalize_weights") &&
        !parse_weight_boolean(config.get(section, "normalize_weights"), "normalize_weights")) {
        throw std::runtime_error("normalize_weights must be true for statistical CLI calculations.");
    }
    if (config.has(section, "zero_first_order_moment")) {
        options.zero_first_order_moment = parse_weight_boolean(
            config.get(section, "zero_first_order_moment"), "zero_first_order_moment");
    }
    return options;
}

inline weighted_moment_core::WeightedMomentOptions parse_uniform_weight_options(const IniConfig& config) {
    return parse_checked_weight_options(config);
}

} // namespace without_sgems

#endif
