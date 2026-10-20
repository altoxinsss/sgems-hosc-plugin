#ifndef WEIGHTED_MOMENT_CORE_H
#define WEIGHTED_MOMENT_CORE_H

#include <cstddef>
#include <string>
#include <vector>

namespace weighted_moment_core {

enum class WeightMode {
    Uniform,
    CoordinateGaussian,
    ValueGaussian,
    SoftHighPorosity,
    AdaptiveUpperTail,
};

struct WeightedMomentOptions {
    WeightMode mode = WeightMode::Uniform;
    double center_x = 0.0;
    double center_y = 0.0;
    double center_z = 0.0;
    // Positive axis half-widths select relative-window coordinate Gaussian
    // weighting. Zero preserves the legacy physical-distance interpretation.
    double coordinate_scale_x = 0.0;
    double coordinate_scale_y = 0.0;
    double coordinate_scale_z = 0.0;
    double bandwidth = 1.0;
    double reference_value = 0.0;
    double min_weight = 0.0;
    bool normalize_weights = true;
    bool zero_first_order_moment = false;

    // Soft one-sided high-value emphasis:
    // (1-strength) + strength*sigmoid((value-threshold)/transition_width).
    double high_porosity_threshold = 0.25;
    double soft_transition_width = 0.02;
    double soft_weight_strength = 1.0;

    // Adaptive upper-tail emphasis.  The support-specific threshold and width
    // are prepared from all values in the current window/template replicate set.
    double adaptive_global_threshold = 0.25;
    double adaptive_mean_offset = 1.0;
    double adaptive_transition_scale = 0.25;
    double adaptive_weight_strength = 1.0;
    double adaptive_std_floor = 1.0e-6;
    double adaptive_context_mean = 0.0;
    double adaptive_context_stddev = 0.0;
    double adaptive_context_threshold = 0.25;
    double adaptive_context_transition_width = 1.0e-6;
    bool adaptive_context_ready = false;

    // False preserves the historical template behaviour in which values[0]
    // supplies the configuration weight.  True uses the geometric mean of
    // the component-level value weights.
    bool joint_geometric_mean = false;
};

struct WeightedSample {
    const double* values = nullptr;
    int value_count = 0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct CumulantPartitionTerm {
    long double coefficient = 0.0;
    std::vector<unsigned int> block_masks;
};

// The combinatorial structure of a joint cumulant depends only on the
// multiplicity vector.  It can therefore be compiled once and reused for
// every local-template replicate set.
struct PreparedCumulantPlan {
    std::vector<int> labels;
    std::vector<unsigned int> required_masks;
    std::vector<CumulantPartitionTerm> terms;
    int order = 0;
};

struct WeightDiagnostics {
    int positive_weight_count = 0;
    double weight_sum = 0.0;
    double weight_square_sum = 0.0;
    double effective_sample_size = 0.0;
};

double compute_weight(const WeightedSample& sample, const WeightedMomentOptions& options);

WeightedMomentOptions prepare_weight_context(
    const WeightedSample* samples,
    int sample_count,
    const WeightedMomentOptions& options);

std::vector<double> compute_univariate_weights(
    const double* values,
    int sample_count,
    const WeightedMomentOptions& options,
    WeightedMomentOptions* resolved_options = nullptr);

const char* weighting_formula_version();
const char* weight_mode_name(WeightMode mode);
std::string weight_scheme_id(const WeightedMomentOptions& options);
std::string weight_options_manifest(const WeightedMomentOptions& options);

PreparedCumulantPlan prepare_cumulant_plan(
    const int* multiplicities,
    int variable_count);

double evaluate_weighted_cumulant(
    const WeightedSample* samples,
    int sample_count,
    const PreparedCumulantPlan& plan,
    const WeightedMomentOptions& options);

WeightDiagnostics compute_weight_diagnostics(
    const WeightedSample* samples,
    int sample_count,
    const WeightedMomentOptions& options);

double weighted_moment(
    const WeightedSample* samples,
    int sample_count,
    const int* multiplicities,
    int variable_count,
    const WeightedMomentOptions& options);

double weighted_cumulant_from_multiplicities(
    const WeightedSample* samples,
    int sample_count,
    const int* multiplicities,
    int variable_count,
    const WeightedMomentOptions& options);

double weighted_univariate_moment(
    const double* values,
    const double* weights,
    int sample_count,
    int order,
    bool normalize_weights = true,
    bool zero_first_order_moment = false);

double weighted_univariate_cumulant(
    const double* values,
    const double* weights,
    int sample_count,
    int order,
    bool normalize_weights = true,
    bool zero_first_order_moment = false);

} // namespace weighted_moment_core

#endif
