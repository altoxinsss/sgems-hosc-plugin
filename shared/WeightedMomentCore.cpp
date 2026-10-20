#include "WeightedMomentCore.h"

#include <cmath>
#include <float.h>
#include <limits>
#include <stdexcept>

namespace weighted_moment_core {
namespace {

const char* const kWeightingFormulaVersion = "2026-10-15.v3-relative-window";

double max_double(double left, double right) {
    return left > right ? left : right;
}

double min_double(double left, double right) {
    return left < right ? left : right;
}

double clamp_unit(double value) {
    return max_double(0.0, min_double(1.0, value));
}

double stable_sigmoid(double value) {
    if(value >= 0.0) {
        return 1.0 / (1.0 + std::exp(-value));
    }
    const double exponent = std::exp(value);
    return exponent / (1.0 + exponent);
}

double binomial(int n, int k) {
    if(k < 0 || k > n) {
        return 0.0;
    }
    if(k == 0 || k == n) {
        return 1.0;
    }
    if(k > n - k) {
        k = n - k;
    }
    long double result = 1.0;
    int i = 0;
    for(i = 1; i <= k; ++i) {
        result *= static_cast<long double>(n - k + i);
        result /= static_cast<long double>(i);
    }
    return static_cast<double>(result);
}

double factorial_value(int n) {
    long double value = 1.0;
    int i = 0;
    for(i = 2; i <= n; ++i) {
        value *= static_cast<long double>(i);
    }
    return static_cast<double>(value);
}

bool nonnegative_finite(double value) {
    return is_finite_number(value) && value >= 0.0;
}

bool value_mode(WeightMode mode) {
    return mode == WeightMode_ValueGaussian || mode == WeightMode_SoftHighPorosity;
}

double component_value_weight(double value, const WeightedMomentOptions& options) {
    const double epsilon = std::numeric_limits<double>::epsilon();
    if(options.mode == WeightMode_ValueGaussian) {
        const double bandwidth = max_double(options.bandwidth, epsilon);
        const double distance = value - options.reference_value;
        return std::exp(-0.5 * distance * distance / (bandwidth * bandwidth));
    }
    if(options.mode == WeightMode_SoftHighPorosity) {
        const double width = max_double(options.soft_transition_width, epsilon);
        const double strength = clamp_unit(options.soft_weight_strength);
        return (1.0 - strength) + strength * stable_sigmoid(
            (value - options.high_porosity_threshold) / width);
    }
    return 1.0;
}

double normalized_axis_distance(double difference, double scale) {
    if(scale <= 0.0) {
        // A singleton axis has zero displacement from its window center.
        return 0.0;
    }
    return difference / scale;
}

void validate_multiplicities(const int* multiplicities, int variable_count) {
    int i = 0;
    if(multiplicities == NULL || variable_count <= 0) {
        throw std::invalid_argument("weighted cumulant requires a non-empty multiplicity index.");
    }
    for(i = 0; i < variable_count; ++i) {
        if(multiplicities[i] < 0) {
            throw std::invalid_argument("weighted cumulant multiplicities must be non-negative.");
        }
    }
}

int joint_order(const int* multiplicities, int variable_count) {
    int result = 0;
    int i = 0;
    for(i = 0; i < variable_count; ++i) {
        result += multiplicities[i];
    }
    return result;
}

std::vector<int> labels_from_multiplicities(const int* multiplicities, int variable_count) {
    std::vector<int> labels;
    int variable = 0;
    labels.reserve(static_cast<size_t>(joint_order(multiplicities, variable_count)));
    for(variable = 0; variable < variable_count; ++variable) {
        int repeat = 0;
        for(repeat = 0; repeat < multiplicities[variable]; ++repeat) {
            labels.push_back(variable);
        }
    }
    return labels;
}

std::vector<double> weighted_variable_centers(
    const WeightedSample* samples,
    int sample_count,
    int variable_count,
    const WeightedMomentOptions& options) {
    std::vector<double> centers(static_cast<size_t>(variable_count), 0.0);
    std::vector<long double> sums(static_cast<size_t>(variable_count), 0.0);
    long double sum_weight = 0.0;
    int sample = 0;
    for(sample = 0; sample < sample_count; ++sample) {
        const double weight = compute_weight(samples[sample], options);
        int variable = 0;
        if(weight <= 0.0) {
            continue;
        }
        if(samples[sample].values == NULL || samples[sample].value_count < variable_count) {
            throw std::invalid_argument("weighted cumulant sample does not contain all requested variables.");
        }
        sum_weight += static_cast<long double>(weight);
        for(variable = 0; variable < variable_count; ++variable) {
            sums[static_cast<size_t>(variable)] += static_cast<long double>(weight) * samples[sample].values[variable];
        }
    }
    if(sum_weight > 0.0) {
        int variable = 0;
        for(variable = 0; variable < variable_count; ++variable) {
            centers[static_cast<size_t>(variable)] = static_cast<double>(
                sums[static_cast<size_t>(variable)] / sum_weight);
        }
    }
    return centers;
}

double product_for_mask(
    const WeightedSample& sample,
    const std::vector<int>& labels,
    int mask,
    const std::vector<double>& centers) {
    long double product = 1.0;
    int i = 0;
    for(i = 0; i < static_cast<int>(labels.size()); ++i) {
        if((mask & (1 << i)) == 0) {
            continue;
        }
        const int variable = labels[static_cast<size_t>(i)];
        const double center = centers.empty() ? 0.0 : centers[static_cast<size_t>(variable)];
        if(sample.values == NULL || variable < 0 || variable >= sample.value_count) {
            throw std::invalid_argument("weighted cumulant sample does not contain all requested variables.");
        }
        product *= static_cast<long double>(sample.values[variable] - center);
    }
    return static_cast<double>(product);
}

void enumerate_partitions_rec(
    int count,
    int next,
    std::vector<int>& block_for_item,
    int block_count,
    std::vector< std::vector<int> >& partitions) {
    int block = 0;
    if(next == count) {
        partitions.push_back(block_for_item);
        return;
    }
    for(block = 0; block < block_count; ++block) {
        block_for_item[static_cast<size_t>(next)] = block;
        enumerate_partitions_rec(count, next + 1, block_for_item, block_count, partitions);
    }
    block_for_item[static_cast<size_t>(next)] = block_count;
    enumerate_partitions_rec(count, next + 1, block_for_item, block_count + 1, partitions);
}

std::vector< std::vector<int> > enumerate_partitions(int count) {
    std::vector< std::vector<int> > partitions;
    if(count > 0) {
        std::vector<int> block_for_item(static_cast<size_t>(count), 0);
        enumerate_partitions_rec(count, 1, block_for_item, 1, partitions);
    }
    return partitions;
}

double weighted_block_moment(
    const WeightedSample* samples,
    int sample_count,
    const std::vector<int>& labels,
    int mask,
    const WeightedMomentOptions& options) {
    long double sum = 0.0;
    long double sum_weight = 0.0;
    int sample = 0;
    int variable_count = 0;
    std::vector<double> centers;
    if(mask == 0) {
        return 1.0;
    }
    for(sample = 0; sample < static_cast<int>(labels.size()); ++sample) {
        if(labels[static_cast<size_t>(sample)] + 1 > variable_count) {
            variable_count = labels[static_cast<size_t>(sample)] + 1;
        }
    }
    if(options.zero_first_order_moment) {
        centers = weighted_variable_centers(samples, sample_count, variable_count, options);
    }
    for(sample = 0; sample < sample_count; ++sample) {
        const double weight = compute_weight(samples[sample], options);
        if(weight <= 0.0) {
            continue;
        }
        sum += static_cast<long double>(weight) * product_for_mask(samples[sample], labels, mask, centers);
        sum_weight += static_cast<long double>(weight);
    }
    if(sum_weight <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return static_cast<double>(options.normalize_weights ? sum / sum_weight : sum);
}

} // namespace

WeightedMomentOptions::WeightedMomentOptions()
    : mode(WeightMode_Uniform),
      center_x(0.0), center_y(0.0), center_z(0.0),
      coordinate_scale_x(0.0), coordinate_scale_y(0.0), coordinate_scale_z(0.0),
      bandwidth(1.0), reference_value(0.0), min_weight(0.0),
      normalize_weights(true), zero_first_order_moment(false),
      high_porosity_threshold(0.25), soft_transition_width(0.02),
      soft_weight_strength(0.8), joint_geometric_mean(false) {
}

WeightedSample::WeightedSample()
    : values(NULL), value_count(0), x(0.0), y(0.0), z(0.0) {
}

WeightDiagnostics::WeightDiagnostics()
    : positive_weight_count(0), weight_sum(0.0), weight_square_sum(0.0), effective_sample_size(0.0) {
}

const char* weighting_formula_version() {
    return kWeightingFormulaVersion;
}

const char* weight_mode_name(WeightMode mode) {
    switch(mode) {
    case WeightMode_Uniform: return "uniform";
    case WeightMode_CoordinateGaussian: return "coordinate_gaussian";
    case WeightMode_SoftHighPorosity: return "soft_high_porosity";
    case WeightMode_ValueGaussian: return "value_gaussian_legacy";
    case WeightMode_AdaptiveUpperTail: return "adaptive_upper_tail_legacy";
    default: return "unknown";
    }
}

bool is_finite_number(double value) {
#if defined(_MSC_VER)
    return _finite(value) != 0;
#else
    return std::isfinite(value);
#endif
}

bool configure_relative_window_coordinate_gaussian(
    WeightedMomentOptions* options,
    double origin_x, double origin_y, double origin_z,
    double spacing_x, double spacing_y, double spacing_z,
    int start_x, int start_y, int start_z,
    int size_x, int size_y, int size_z) {
    if(options == NULL || size_x <= 0 || size_y <= 0 || size_z <= 0 ||
       !is_finite_number(origin_x) || !is_finite_number(origin_y) || !is_finite_number(origin_z) ||
       !is_finite_number(spacing_x) || !is_finite_number(spacing_y) || !is_finite_number(spacing_z)) {
        return false;
    }
    options->center_x = origin_x + (static_cast<double>(start_x) +
        0.5 * static_cast<double>(size_x - 1)) * spacing_x;
    options->center_y = origin_y + (static_cast<double>(start_y) +
        0.5 * static_cast<double>(size_y - 1)) * spacing_y;
    options->center_z = origin_z + (static_cast<double>(start_z) +
        0.5 * static_cast<double>(size_z - 1)) * spacing_z;
    options->coordinate_scale_x = 0.5 * static_cast<double>(size_x - 1) * std::fabs(spacing_x);
    options->coordinate_scale_y = 0.5 * static_cast<double>(size_y - 1) * std::fabs(spacing_y);
    options->coordinate_scale_z = 0.5 * static_cast<double>(size_z - 1) * std::fabs(spacing_z);
    return true;
}

double compute_weight(const WeightedSample& sample, const WeightedMomentOptions& options) {
    const double epsilon = std::numeric_limits<double>::epsilon();
    double weight = 1.0;
    if(options.mode == WeightMode_CoordinateGaussian) {
        const double dx = normalized_axis_distance(sample.x - options.center_x, options.coordinate_scale_x);
        const double dy = normalized_axis_distance(sample.y - options.center_y, options.coordinate_scale_y);
        const double dz = normalized_axis_distance(sample.z - options.center_z, options.coordinate_scale_z);
        const double bandwidth = max_double(options.bandwidth, epsilon);
        weight = std::exp(-0.5 * (dx * dx + dy * dy + dz * dz) / (bandwidth * bandwidth));
    }
    else if(value_mode(options.mode)) {
        if(sample.values == NULL || sample.value_count <= 0) {
            return 0.0;
        }
        if(options.joint_geometric_mean && sample.value_count > 1) {
            long double log_sum = 0.0;
            int i = 0;
            for(i = 0; i < sample.value_count; ++i) {
                const double component = max_double(component_value_weight(sample.values[i], options),
                    std::numeric_limits<double>::min());
                log_sum += std::log(component);
            }
            weight = std::exp(static_cast<double>(log_sum / sample.value_count));
        }
        else {
            weight = component_value_weight(sample.values[0], options);
        }
    }
    if(!nonnegative_finite(weight) || weight < options.min_weight) {
        return 0.0;
    }
    return weight;
}

std::vector<double> compute_univariate_weights(
    const double* values, int sample_count, const double* x, const double* y, const double* z,
    const WeightedMomentOptions& options) {
    std::vector<double> weights;
    int i = 0;
    if(values == NULL || sample_count <= 0) {
        throw std::invalid_argument("compute_univariate_weights received invalid input.");
    }
    weights.reserve(static_cast<size_t>(sample_count));
    for(i = 0; i < sample_count; ++i) {
        WeightedSample sample;
        sample.values = values + i;
        sample.value_count = 1;
        sample.x = x == NULL ? 0.0 : x[i];
        sample.y = y == NULL ? 0.0 : y[i];
        sample.z = z == NULL ? 0.0 : z[i];
        weights.push_back(compute_weight(sample, options));
    }
    return weights;
}

WeightDiagnostics compute_weight_diagnostics(const double* weights, int sample_count) {
    WeightDiagnostics result;
    long double sum = 0.0;
    long double squared_sum = 0.0;
    int i = 0;
    if(weights == NULL || sample_count <= 0) {
        return result;
    }
    for(i = 0; i < sample_count; ++i) {
        const double weight = weights[i];
        if(!nonnegative_finite(weight) || weight <= 0.0) {
            continue;
        }
        ++result.positive_weight_count;
        sum += weight;
        squared_sum += static_cast<long double>(weight) * weight;
    }
    result.weight_sum = static_cast<double>(sum);
    result.weight_square_sum = static_cast<double>(squared_sum);
    if(squared_sum > 0.0) {
        result.effective_sample_size = static_cast<double>(sum * sum / squared_sum);
    }
    return result;
}

double weighted_moment(
    const WeightedSample* samples, int sample_count, const int* multiplicities,
    int variable_count, const WeightedMomentOptions& options) {
    std::vector<int> labels;
    int order = 0;
    validate_multiplicities(multiplicities, variable_count);
    labels = labels_from_multiplicities(multiplicities, variable_count);
    order = static_cast<int>(labels.size());
    if(samples == NULL || sample_count <= 0) {
        throw std::invalid_argument("weighted moment requires at least one sample.");
    }
    if(order == 0) {
        return 1.0;
    }
    if(order >= static_cast<int>(sizeof(int) * 8 - 1)) {
        throw std::invalid_argument("weighted moment order is too large for masks.");
    }
    return weighted_block_moment(samples, sample_count, labels, (1 << order) - 1, options);
}

double weighted_cumulant_from_multiplicities(
    const WeightedSample* samples, int sample_count, const int* multiplicities,
    int variable_count, const WeightedMomentOptions& options) {
    std::vector<int> labels;
    int order = 0;
    std::vector< std::vector<int> > partitions;
    long double cumulant = 0.0;
    size_t partition_index = 0;
    if(samples == NULL || sample_count <= 0) {
        throw std::invalid_argument("weighted cumulant requires at least one sample.");
    }
    validate_multiplicities(multiplicities, variable_count);
    labels = labels_from_multiplicities(multiplicities, variable_count);
    order = static_cast<int>(labels.size());
    if(order == 0) {
        return 0.0;
    }
    if(order >= static_cast<int>(sizeof(int) * 8 - 1)) {
        throw std::invalid_argument("weighted cumulant order is too large for partition masks.");
    }
    partitions = enumerate_partitions(order);
    for(partition_index = 0; partition_index < partitions.size(); ++partition_index) {
        const std::vector<int>& partition = partitions[partition_index];
        int block_count = 0;
        int item = 0;
        long double term;
        for(item = 0; item < order; ++item) {
            if(partition[static_cast<size_t>(item)] + 1 > block_count) {
                block_count = partition[static_cast<size_t>(item)] + 1;
            }
        }
        term = factorial_value(block_count - 1);
        if((block_count - 1) % 2 != 0) {
            term = -term;
        }
        for(item = 0; item < block_count; ++item) {
            int mask = 0;
            int index = 0;
            for(index = 0; index < order; ++index) {
                if(partition[static_cast<size_t>(index)] == item) {
                    mask |= (1 << index);
                }
            }
            term *= weighted_block_moment(samples, sample_count, labels, mask, options);
        }
        cumulant += term;
    }
    return static_cast<double>(cumulant);
}

double weighted_univariate_moment(
    const double* values, const double* weights, int sample_count, int order,
    bool normalize_weights, bool zero_first_order_moment) {
    long double weighted_sum = 0.0;
    long double sum_weight = 0.0;
    long double weighted_mean = 0.0;
    double center = 0.0;
    int i = 0;
    if(values == NULL || sample_count <= 0 || order < 0) {
        throw std::invalid_argument("weighted univariate moment received invalid input.");
    }
    if(order == 0) {
        return 1.0;
    }
    if(zero_first_order_moment) {
        for(i = 0; i < sample_count; ++i) {
            const double weight = weights == NULL ? 1.0 : weights[i];
            if(!nonnegative_finite(weight) || weight <= 0.0 || !is_finite_number(values[i])) {
                continue;
            }
            weighted_mean += static_cast<long double>(weight) * values[i];
            sum_weight += weight;
        }
        if(sum_weight <= 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        center = static_cast<double>(weighted_mean / sum_weight);
    }
    sum_weight = 0.0;
    for(i = 0; i < sample_count; ++i) {
        const double weight = weights == NULL ? 1.0 : weights[i];
        if(!nonnegative_finite(weight) || weight <= 0.0 || !is_finite_number(values[i])) {
            continue;
        }
        weighted_sum += static_cast<long double>(weight) * std::pow(values[i] - center, order);
        sum_weight += weight;
    }
    if(sum_weight <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return static_cast<double>(normalize_weights ? weighted_sum / sum_weight : weighted_sum);
}

double weighted_univariate_cumulant(
    const double* values, const double* weights, int sample_count, int order,
    bool normalize_weights, bool zero_first_order_moment) {
    std::vector<double> moments;
    std::vector<double> cumulants;
    int n = 0;
    if(order <= 0) {
        return 0.0;
    }
    moments.assign(static_cast<size_t>(order) + 1, 0.0);
    cumulants.assign(static_cast<size_t>(order) + 1, 0.0);
    for(n = 1; n <= order; ++n) {
        moments[static_cast<size_t>(n)] = weighted_univariate_moment(
            values, weights, sample_count, n, normalize_weights, zero_first_order_moment);
    }
    cumulants[1] = moments[1];
    for(n = 2; n <= order; ++n) {
        long double value = moments[static_cast<size_t>(n)];
        int j = 0;
        for(j = 1; j < n; ++j) {
            value -= static_cast<long double>(binomial(n - 1, j - 1)) *
                cumulants[static_cast<size_t>(j)] * moments[static_cast<size_t>(n - j)];
        }
        cumulants[static_cast<size_t>(n)] = static_cast<double>(value);
    }
    return cumulants[static_cast<size_t>(order)];
}

} // namespace weighted_moment_core
