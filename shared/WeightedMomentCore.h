#ifndef WITHSGEMS_SHARED_WEIGHTED_MOMENT_CORE_H
#define WITHSGEMS_SHARED_WEIGHTED_MOMENT_CORE_H

#include <vector>

namespace weighted_moment_core {

// The two legacy modes remain representable only so old parameter files can
// be rejected or described explicitly by the plugin contract.  New UI paths
// expose Uniform, CoordinateGaussian, and SoftHighPorosity only.
enum WeightMode {
    WeightMode_Uniform,
    WeightMode_CoordinateGaussian,
    WeightMode_ValueGaussian,
    WeightMode_SoftHighPorosity,
    WeightMode_AdaptiveUpperTail
};

struct WeightedMomentOptions {
    WeightMode mode;
    double center_x;
    double center_y;
    double center_z;
    // Positive scales enable the formal relative-window Gaussian.  A zero
    // scale on a singleton axis means that axis contributes zero distance.
    double coordinate_scale_x;
    double coordinate_scale_y;
    double coordinate_scale_z;
    double bandwidth;
    double reference_value;
    double min_weight;
    bool normalize_weights;
    bool zero_first_order_moment;
    double high_porosity_threshold;
    double soft_transition_width;
    double soft_weight_strength;
    bool joint_geometric_mean;

    WeightedMomentOptions();
};

struct WeightedSample {
    const double* values;
    int value_count;
    double x;
    double y;
    double z;

    WeightedSample();
};

struct WeightDiagnostics {
    int positive_weight_count;
    double weight_sum;
    double weight_square_sum;
    double effective_sample_size;

    WeightDiagnostics();
};

const char* weighting_formula_version();
const char* weight_mode_name(WeightMode mode);

bool is_finite_number(double value);
// Configures the relative-window coordinate Gaussian used by HOSC_WINDOW.
// The center follows the actual source-grid geometry; singleton axes receive
// a zero scale so they make no contribution to the normalized distance.
bool configure_relative_window_coordinate_gaussian(
    WeightedMomentOptions* options,
    double origin_x, double origin_y, double origin_z,
    double spacing_x, double spacing_y, double spacing_z,
    int start_x, int start_y, int start_z,
    int size_x, int size_y, int size_z);
double compute_weight(const WeightedSample& sample, const WeightedMomentOptions& options);

// Computes non-normalized contributions. Callers can pass the same vector to
// every requested moment order; weighted_* functions normalize internally
// when normalize_weights is true.
std::vector<double> compute_univariate_weights(
    const double* values,
    int sample_count,
    const double* x,
    const double* y,
    const double* z,
    const WeightedMomentOptions& options);

WeightDiagnostics compute_weight_diagnostics(
    const double* weights,
    int sample_count);

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
    bool normalize_weights,
    bool zero_first_order_moment);

double weighted_univariate_cumulant(
    const double* values,
    const double* weights,
    int sample_count,
    int order,
    bool normalize_weights,
    bool zero_first_order_moment);

} // namespace weighted_moment_core

#endif
