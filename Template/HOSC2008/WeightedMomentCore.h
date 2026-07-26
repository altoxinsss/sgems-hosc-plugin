#ifndef WEIGHTED_MOMENT_CORE_H
#define WEIGHTED_MOMENT_CORE_H

namespace weighted_moment_core {

enum WeightMode {
    WeightMode_Uniform,
    WeightMode_CoordinateGaussian,
    WeightMode_ValueGaussian
};

struct WeightedMomentOptions {
    WeightMode mode;
    double center_x;
    double center_y;
    double center_z;
    double bandwidth;
    double reference_value;
    double min_weight;
    bool normalize_weights;
    bool zero_first_order_moment;

    WeightedMomentOptions()
        : mode(WeightMode_Uniform),
          center_x(0.0),
          center_y(0.0),
          center_z(0.0),
          bandwidth(1.0),
          reference_value(0.0),
          min_weight(0.0),
          normalize_weights(true),
          zero_first_order_moment(false) {}
};

struct WeightedSample {
    const double* values;
    int value_count;
    double x;
    double y;
    double z;

    WeightedSample()
        : values(0),
          value_count(0),
          x(0.0),
          y(0.0),
          z(0.0) {}
};

double compute_weight(const WeightedSample& sample, const WeightedMomentOptions& options);

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

}

#endif
