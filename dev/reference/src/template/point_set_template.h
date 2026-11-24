#ifndef WITHOUT_SGEMS_POINT_SET_TEMPLATE_H
#define WITHOUT_SGEMS_POINT_SET_TEMPLATE_H

#include "config.h"
#include "WeightedMomentCore.h"
#include <filesystem>
#include <vector>

namespace without_sgems {
// Entry point separate from the integer/dense-grid compatibility path.
int run_point_set_template(const IniConfig& config,
    const std::filesystem::path& config_path);

// Shared statistical evaluation only: sampling remains owned by each entry
// point. Values are raw on entry and may be centered for a cumulant target.
double evaluate_complete_template_samples(std::vector<std::vector<double>>& values,
    const std::vector<double>& tuple_weights, const std::vector<int>& nod,
    bool moment, bool bell, const weighted_moment_core::WeightedMomentOptions& options);
}

#endif
