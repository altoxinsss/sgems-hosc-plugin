#include "WeightedMomentCore.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace weighted_moment_core {
namespace {

// v2 records the formal relative-window coordinate Gaussian definition.
// Coordinate offsets are divided by each axis half-width before sigma_s
// (stored in bandwidth) is applied; this is intentionally distinct from the
// legacy physical-distance coordinate Gaussian recorded under v1.
const char* const kWeightingFormulaVersion = "2026-10-10.v2";

double clamp_unit(double value) {
    return std::max(0.0, std::min(1.0, value));
}

double stable_sigmoid(double value) {
    if (value >= 0.0) {
        return 1.0 / (1.0 + std::exp(-value));
    }
    const double exponential = std::exp(value);
    return exponential / (1.0 + exponential);
}

bool is_value_based_mode(WeightMode mode) {
    return mode == WeightMode::ValueGaussian ||
        mode == WeightMode::SoftHighPorosity ||
        mode == WeightMode::AdaptiveUpperTail;
}

double component_value_weight(double value, const WeightedMomentOptions& options) {
    if (options.mode == WeightMode::ValueGaussian) {
        const double bandwidth = std::max(options.bandwidth, std::numeric_limits<double>::epsilon());
        const double difference = value - options.reference_value;
        return std::exp(-0.5 * difference * difference / (bandwidth * bandwidth));
    }
    if (options.mode == WeightMode::SoftHighPorosity) {
        const double width = std::max(options.soft_transition_width, std::numeric_limits<double>::epsilon());
        const double strength = clamp_unit(options.soft_weight_strength);
        const double membership = stable_sigmoid((value - options.high_porosity_threshold) / width);
        return (1.0 - strength) + strength * membership;
    }
    if (options.mode == WeightMode::AdaptiveUpperTail) {
        if (!options.adaptive_context_ready) {
            throw std::invalid_argument(
                "adaptive_upper_tail requires a prepared support context; "
                "use prepare_weight_context or compute_univariate_weights.");
        }
        const double width = std::max(
            options.adaptive_context_transition_width,
            std::numeric_limits<double>::epsilon());
        const double strength = clamp_unit(options.adaptive_weight_strength);
        const double membership = stable_sigmoid(
            (value - options.adaptive_context_threshold) / width);
        return (1.0 - strength) + strength * membership;
    }
    return 1.0;
}

double binomial(int n, int k) {
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

void validate_multiplicities(const int* multiplicities, int variable_count) {
    if (multiplicities == nullptr || variable_count <= 0) {
        throw std::invalid_argument("weighted cumulant requires a non-empty multiplicity index.");
    }
    for (int i = 0; i < variable_count; ++i) {
        if (multiplicities[i] < 0) {
            throw std::invalid_argument("weighted cumulant multiplicities must be non-negative.");
        }
    }
}

int joint_order(const int* multiplicities, int variable_count) {
    int order = 0;
    for (int i = 0; i < variable_count; ++i) {
        order += multiplicities[i];
    }
    return order;
}

std::vector<double> weighted_variable_centers(
    const WeightedSample* samples,
    int sample_count,
    int variable_count,
    const WeightedMomentOptions& options) {
    std::vector<double> centers(static_cast<size_t>(variable_count), 0.0);
    std::vector<long double> sums(static_cast<size_t>(variable_count), 0.0);
    long double weight_sum = 0.0;

    for (int i = 0; i < sample_count; ++i) {
        const double weight = compute_weight(samples[i], options);
        if (weight <= 0.0) {
            continue;
        }
        if (samples[i].values == nullptr || samples[i].value_count < variable_count) {
            throw std::invalid_argument("weighted cumulant sample does not contain all requested variables.");
        }
        weight_sum += static_cast<long double>(weight);
        for (int variable = 0; variable < variable_count; ++variable) {
            sums[static_cast<size_t>(variable)] +=
                static_cast<long double>(weight) * samples[i].values[variable];
        }
    }

    if (weight_sum <= 0.0) {
        return centers;
    }
    for (int variable = 0; variable < variable_count; ++variable) {
        centers[static_cast<size_t>(variable)] =
            static_cast<double>(sums[static_cast<size_t>(variable)] / weight_sum);
    }
    return centers;
}

double product_for_mask(
    const WeightedSample& sample,
    const std::vector<int>& labels,
    int mask,
    const std::vector<double>& centers) {
    long double product = 1.0;
    for (int i = 0; i < static_cast<int>(labels.size()); ++i) {
        if ((mask & (1 << i)) == 0) {
            continue;
        }
        const int variable = labels[static_cast<size_t>(i)];
        if (sample.values == nullptr || variable < 0 || variable >= sample.value_count) {
            throw std::invalid_argument("weighted cumulant sample does not contain all requested variables.");
        }
        const double center = centers.empty() ? 0.0 : centers[static_cast<size_t>(variable)];
        product *= static_cast<long double>(sample.values[variable] - center);
    }
    return static_cast<double>(product);
}

std::vector<int> labels_from_multiplicities(const int* multiplicities, int variable_count) {
    std::vector<int> labels;
    labels.reserve(static_cast<size_t>(joint_order(multiplicities, variable_count)));
    for (int variable = 0; variable < variable_count; ++variable) {
        for (int repeat = 0; repeat < multiplicities[variable]; ++repeat) {
            labels.push_back(variable);
        }
    }
    return labels;
}

void enumerate_partitions_rec(
    int n,
    int next,
    std::vector<int>& block_for_item,
    int block_count,
    std::vector<std::vector<int>>& partitions) {
    if (next == n) {
        partitions.push_back(block_for_item);
        return;
    }
    for (int b = 0; b < block_count; ++b) {
        block_for_item[static_cast<size_t>(next)] = b;
        enumerate_partitions_rec(n, next + 1, block_for_item, block_count, partitions);
    }
    block_for_item[static_cast<size_t>(next)] = block_count;
    enumerate_partitions_rec(n, next + 1, block_for_item, block_count + 1, partitions);
}

std::vector<std::vector<int>> enumerate_partitions(int n) {
    std::vector<std::vector<int>> partitions;
    if (n <= 0) {
        return partitions;
    }
    std::vector<int> block_for_item(static_cast<size_t>(n), 0);
    enumerate_partitions_rec(n, 1, block_for_item, 1, partitions);
    return partitions;
}

double weighted_block_moment(
    const WeightedSample* samples,
    int sample_count,
    const std::vector<int>& labels,
    int mask,
    const WeightedMomentOptions& options) {
    if (mask == 0) {
        return 1.0;
    }
    int variable_count = 0;
    for (size_t i = 0; i < labels.size(); ++i) {
        variable_count = std::max(variable_count, labels[i] + 1);
    }
    const std::vector<double> centers = options.zero_first_order_moment
        ? weighted_variable_centers(samples, sample_count, variable_count, options)
        : std::vector<double>();
    long double weighted_sum = 0.0;
    long double weight_sum = 0.0;
    for (int i = 0; i < sample_count; ++i) {
        const double w = compute_weight(samples[i], options);
        if (w <= 0.0) {
            continue;
        }
        weighted_sum += static_cast<long double>(w) * product_for_mask(samples[i], labels, mask, centers);
        weight_sum += static_cast<long double>(w);
    }
    if (weight_sum <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return static_cast<double>(options.normalize_weights ? weighted_sum / weight_sum : weighted_sum);
}

double weighted_block_moment_with_centers(
    const WeightedSample* samples,
    int sample_count,
    const std::vector<int>& labels,
    unsigned int mask,
    const std::vector<double>& centers,
    const WeightedMomentOptions& options) {
    if (mask == 0U) {
        return 1.0;
    }
    long double weighted_sum = 0.0;
    long double weight_sum = 0.0;
    for (int i = 0; i < sample_count; ++i) {
        const double weight = compute_weight(samples[i], options);
        if (weight <= 0.0) {
            continue;
        }
        weighted_sum += static_cast<long double>(weight) *
            product_for_mask(samples[i], labels, static_cast<int>(mask), centers);
        weight_sum += static_cast<long double>(weight);
    }
    if (weight_sum <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return static_cast<double>(options.normalize_weights ? weighted_sum / weight_sum : weighted_sum);
}

} // namespace

double compute_weight(const WeightedSample& sample, const WeightedMomentOptions& options) {
    double weight = 1.0;
    const double bandwidth = std::max(options.bandwidth, std::numeric_limits<double>::epsilon());
    if (options.mode == WeightMode::CoordinateGaussian) {
        const double dx = sample.x - options.center_x;
        const double dy = sample.y - options.center_y;
        const double dz = sample.z - options.center_z;
        const bool relative_coordinates =
            options.coordinate_scale_x > 0.0 && options.coordinate_scale_y > 0.0 && options.coordinate_scale_z > 0.0;
        const double sx = relative_coordinates ? options.coordinate_scale_x : 1.0;
        const double sy = relative_coordinates ? options.coordinate_scale_y : 1.0;
        const double sz = relative_coordinates ? options.coordinate_scale_z : 1.0;
        const double d2 = dx * dx / (sx * sx) + dy * dy / (sy * sy) + dz * dz / (sz * sz);
        weight = std::exp(-0.5 * d2 / (bandwidth * bandwidth));
    }
    else if (is_value_based_mode(options.mode)) {
        if (sample.values == nullptr || sample.value_count <= 0) {
            return 0.0;
        }
        if (options.joint_geometric_mean && sample.value_count > 1) {
            long double log_sum = 0.0;
            for (int i = 0; i < sample.value_count; ++i) {
                const double component_weight = std::max(
                    component_value_weight(sample.values[i], options),
                    std::numeric_limits<double>::min());
                log_sum += std::log(component_weight);
            }
            weight = std::exp(static_cast<double>(
                log_sum / static_cast<long double>(sample.value_count)));
        }
        else {
            weight = component_value_weight(sample.values[0], options);
        }
    }
    if (weight < options.min_weight) {
        weight = 0.0;
    }
    return weight;
}

WeightedMomentOptions prepare_weight_context(
    const WeightedSample* samples,
    int sample_count,
    const WeightedMomentOptions& options) {
    WeightedMomentOptions resolved = options;
    if (options.mode != WeightMode::AdaptiveUpperTail) {
        return resolved;
    }
    if (samples == nullptr || sample_count <= 0) {
        throw std::invalid_argument("adaptive_upper_tail requires at least one sample.");
    }

    long double sum = 0.0;
    long double sum_squared = 0.0;
    long long count = 0;
    for (int i = 0; i < sample_count; ++i) {
        if (samples[i].values == nullptr || samples[i].value_count <= 0) {
            continue;
        }
        for (int j = 0; j < samples[i].value_count; ++j) {
            const double value = samples[i].values[j];
            if (!std::isfinite(value)) {
                continue;
            }
            sum += static_cast<long double>(value);
            sum_squared += static_cast<long double>(value) * static_cast<long double>(value);
            ++count;
        }
    }
    if (count <= 0) {
        throw std::invalid_argument("adaptive_upper_tail support contains no finite values.");
    }

    const long double mean = sum / static_cast<long double>(count);
    const long double variance = std::max(
        static_cast<long double>(0.0),
        sum_squared / static_cast<long double>(count) - mean * mean);
    const double standard_deviation = std::sqrt(static_cast<double>(variance));
    const double standard_deviation_for_width = std::max(
        standard_deviation,
        std::max(options.adaptive_std_floor, std::numeric_limits<double>::epsilon()));

    resolved.adaptive_context_mean = static_cast<double>(mean);
    resolved.adaptive_context_stddev = standard_deviation;
    resolved.adaptive_context_threshold = std::max(
        options.adaptive_global_threshold,
        resolved.adaptive_context_mean + options.adaptive_mean_offset * standard_deviation);
    resolved.adaptive_context_transition_width = std::max(
        options.adaptive_transition_scale * standard_deviation_for_width,
        std::max(options.adaptive_std_floor, std::numeric_limits<double>::epsilon()));
    resolved.adaptive_context_ready = true;
    return resolved;
}

std::vector<double> compute_univariate_weights(
    const double* values,
    int sample_count,
    const WeightedMomentOptions& options,
    WeightedMomentOptions* resolved_options) {
    if (values == nullptr || sample_count <= 0) {
        throw std::invalid_argument("compute_univariate_weights received invalid input.");
    }
    std::vector<WeightedSample> samples(static_cast<size_t>(sample_count));
    for (int i = 0; i < sample_count; ++i) {
        samples[static_cast<size_t>(i)].values = values + i;
        samples[static_cast<size_t>(i)].value_count = 1;
    }
    const WeightedMomentOptions resolved = prepare_weight_context(samples.data(), sample_count, options);
    if (resolved_options != nullptr) {
        *resolved_options = resolved;
    }
    std::vector<double> weights(static_cast<size_t>(sample_count), 0.0);
    for (int i = 0; i < sample_count; ++i) {
        weights[static_cast<size_t>(i)] = compute_weight(samples[static_cast<size_t>(i)], resolved);
    }
    return weights;
}

const char* weighting_formula_version() {
    return kWeightingFormulaVersion;
}

const char* weight_mode_name(WeightMode mode) {
    switch (mode) {
    case WeightMode::CoordinateGaussian:
        return "coordinate_gaussian";
    case WeightMode::ValueGaussian:
        return "fixed_value_gaussian";
    case WeightMode::SoftHighPorosity:
        return "soft_high_porosity";
    case WeightMode::AdaptiveUpperTail:
        return "adaptive_upper_tail";
    case WeightMode::Uniform:
    default:
        return "uniform";
    }
}

std::string weight_scheme_id(const WeightedMomentOptions& options) {
    std::string id = std::string(weight_mode_name(options.mode)) + ".v1";
    if (options.joint_geometric_mean && is_value_based_mode(options.mode)) {
        id += "+joint_geometric_mean.v1";
    }
    return id;
}

std::string weight_options_manifest(const WeightedMomentOptions& options) {
    std::ostringstream stream;
    stream.precision(17);
    stream << "weighting_formula_version=" << weighting_formula_version() << "\n";
    stream << "weighting_scheme_id=" << weight_scheme_id(options) << "\n";
    stream << "weight_mode=" << weight_mode_name(options.mode) << "\n";
    stream << "enable_coordinate_gaussian="
           << (options.mode == WeightMode::CoordinateGaussian ? "true" : "false") << "\n";
    stream << "enable_fixed_value_gaussian="
           << (options.mode == WeightMode::ValueGaussian ? "true" : "false") << "\n";
    stream << "enable_soft_high_porosity="
           << (options.mode == WeightMode::SoftHighPorosity ? "true" : "false") << "\n";
    stream << "enable_adaptive_upper_tail="
           << (options.mode == WeightMode::AdaptiveUpperTail ? "true" : "false") << "\n";
    stream << "enable_joint_geometric_mean=" << (options.joint_geometric_mean ? "true" : "false") << "\n";
    stream << "normalize_weights=" << (options.normalize_weights ? "true" : "false") << "\n";
    stream << "zero_first_order_moment=" << (options.zero_first_order_moment ? "true" : "false") << "\n";
    stream << "min_weight=" << options.min_weight << "\n";
    stream << "reference_value=" << options.reference_value << "\n";
    stream << "bandwidth=" << options.bandwidth << "\n";
    stream << "coordinate_scale_x=" << options.coordinate_scale_x << "\n";
    stream << "coordinate_scale_y=" << options.coordinate_scale_y << "\n";
    stream << "coordinate_scale_z=" << options.coordinate_scale_z << "\n";
    stream << "high_porosity_threshold=" << options.high_porosity_threshold << "\n";
    stream << "soft_transition_width=" << options.soft_transition_width << "\n";
    stream << "soft_weight_strength=" << options.soft_weight_strength << "\n";
    stream << "adaptive_global_threshold=" << options.adaptive_global_threshold << "\n";
    stream << "adaptive_mean_offset=" << options.adaptive_mean_offset << "\n";
    stream << "adaptive_transition_scale=" << options.adaptive_transition_scale << "\n";
    stream << "adaptive_weight_strength=" << options.adaptive_weight_strength << "\n";
    stream << "adaptive_std_floor=" << options.adaptive_std_floor << "\n";
    return stream.str();
}

PreparedCumulantPlan prepare_cumulant_plan(
    const int* multiplicities,
    int variable_count) {
    validate_multiplicities(multiplicities, variable_count);
    PreparedCumulantPlan plan;
    plan.labels = labels_from_multiplicities(multiplicities, variable_count);
    plan.order = static_cast<int>(plan.labels.size());
    if (plan.order == 0) {
        return plan;
    }
    if (plan.order >= static_cast<int>(sizeof(unsigned int) * 8U - 1U)) {
        throw std::invalid_argument("weighted cumulant order is too large for partition masks.");
    }

    const std::vector<std::vector<int>> partitions = enumerate_partitions(plan.order);
    for (size_t p = 0; p < partitions.size(); ++p) {
        const std::vector<int>& partition = partitions[p];
        const int block_count = *std::max_element(partition.begin(), partition.end()) + 1;
        CumulantPartitionTerm term;
        term.coefficient = std::tgamma(static_cast<long double>(block_count));
        if ((block_count - 1) % 2 != 0) {
            term.coefficient = -term.coefficient;
        }
        term.block_masks.reserve(static_cast<size_t>(block_count));
        for (int block = 0; block < block_count; ++block) {
            unsigned int mask = 0U;
            for (int item = 0; item < plan.order; ++item) {
                if (partition[static_cast<size_t>(item)] == block) {
                    mask |= (1U << item);
                }
            }
            term.block_masks.push_back(mask);
            if (std::find(plan.required_masks.begin(), plan.required_masks.end(), mask) ==
                plan.required_masks.end()) {
                plan.required_masks.push_back(mask);
            }
        }
        plan.terms.push_back(term);
    }
    std::sort(plan.required_masks.begin(), plan.required_masks.end());
    return plan;
}

double evaluate_weighted_cumulant(
    const WeightedSample* samples,
    int sample_count,
    const PreparedCumulantPlan& plan,
    const WeightedMomentOptions& options) {
    if (samples == nullptr || sample_count <= 0) {
        throw std::invalid_argument("weighted cumulant requires at least one sample.");
    }
    if (plan.order == 0) {
        return 0.0;
    }
    if (static_cast<int>(plan.labels.size()) != plan.order || plan.terms.empty()) {
        throw std::invalid_argument("Prepared cumulant plan is incomplete.");
    }
    const WeightedMomentOptions resolved_options = prepare_weight_context(samples, sample_count, options);
    int variable_count = 0;
    for (size_t i = 0; i < plan.labels.size(); ++i) {
        variable_count = std::max(variable_count, plan.labels[i] + 1);
    }
    const std::vector<double> centers = resolved_options.zero_first_order_moment
        ? weighted_variable_centers(samples, sample_count, variable_count, resolved_options)
        : std::vector<double>();

    std::map<unsigned int, double> moments;
    for (size_t i = 0; i < plan.required_masks.size(); ++i) {
        const unsigned int mask = plan.required_masks[i];
        const double moment = weighted_block_moment_with_centers(
            samples, sample_count, plan.labels, mask, centers, resolved_options);
        if (!std::isfinite(moment)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        moments[mask] = moment;
    }

    long double cumulant = 0.0;
    for (size_t i = 0; i < plan.terms.size(); ++i) {
        long double term = plan.terms[i].coefficient;
        for (size_t j = 0; j < plan.terms[i].block_masks.size(); ++j) {
            term *= moments[plan.terms[i].block_masks[j]];
        }
        cumulant += term;
    }
    return static_cast<double>(cumulant);
}

WeightDiagnostics compute_weight_diagnostics(
    const WeightedSample* samples,
    int sample_count,
    const WeightedMomentOptions& options) {
    if (samples == nullptr || sample_count <= 0) {
        throw std::invalid_argument("weight diagnostics requires at least one sample.");
    }
    const WeightedMomentOptions resolved_options = prepare_weight_context(samples, sample_count, options);
    long double sum = 0.0;
    long double sum_squares = 0.0;
    WeightDiagnostics diagnostics;
    for (int i = 0; i < sample_count; ++i) {
        const double weight = compute_weight(samples[i], resolved_options);
        if (weight <= 0.0 || !std::isfinite(weight)) {
            continue;
        }
        ++diagnostics.positive_weight_count;
        sum += static_cast<long double>(weight);
        sum_squares += static_cast<long double>(weight) * weight;
    }
    diagnostics.weight_sum = static_cast<double>(sum);
    diagnostics.weight_square_sum = static_cast<double>(sum_squares);
    if (sum_squares > 0.0) {
        diagnostics.effective_sample_size = static_cast<double>(sum * sum / sum_squares);
    }
    return diagnostics;
}

double weighted_moment(
    const WeightedSample* samples,
    int sample_count,
    const int* multiplicities,
    int variable_count,
    const WeightedMomentOptions& options) {
    if (samples == nullptr || sample_count <= 0) {
        throw std::invalid_argument("weighted moment requires at least one sample.");
    }
    validate_multiplicities(multiplicities, variable_count);
    const WeightedMomentOptions resolved_options = prepare_weight_context(samples, sample_count, options);
    const std::vector<int> labels = labels_from_multiplicities(multiplicities, variable_count);
    if (labels.empty()) {
        return 1.0;
    }
    const int full_mask = (1 << static_cast<int>(labels.size())) - 1;
    return weighted_block_moment(samples, sample_count, labels, full_mask, resolved_options);
}

double weighted_cumulant_from_multiplicities(
    const WeightedSample* samples,
    int sample_count,
    const int* multiplicities,
    int variable_count,
    const WeightedMomentOptions& options) {
    if (samples == nullptr || sample_count <= 0) {
        throw std::invalid_argument("weighted cumulant requires at least one sample.");
    }
    const PreparedCumulantPlan plan = prepare_cumulant_plan(multiplicities, variable_count);
    return evaluate_weighted_cumulant(samples, sample_count, plan, options);
}

double weighted_univariate_moment(
    const double* values,
    const double* weights,
    int sample_count,
    int order,
    bool normalize_weights,
    bool zero_first_order_moment) {
    if (values == nullptr || sample_count <= 0 || order < 0) {
        throw std::invalid_argument("weighted univariate moment received invalid input.");
    }
    if (order == 0) {
        return 1.0;
    }
    long double weighted_sum = 0.0;
    long double weight_sum = 0.0;
    long double weighted_mean_sum = 0.0;
    if (zero_first_order_moment) {
        for (int i = 0; i < sample_count; ++i) {
            const double w = (weights == nullptr) ? 1.0 : weights[i];
            if (w <= 0.0) {
                continue;
            }
            weighted_mean_sum += static_cast<long double>(w) * values[i];
            weight_sum += static_cast<long double>(w);
        }
        if (weight_sum <= 0.0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    const double center = zero_first_order_moment
        ? static_cast<double>(weighted_mean_sum / weight_sum)
        : 0.0;
    weight_sum = 0.0;
    for (int i = 0; i < sample_count; ++i) {
        const double w = (weights == nullptr) ? 1.0 : weights[i];
        if (w <= 0.0) {
            continue;
        }
        weighted_sum += static_cast<long double>(w) * std::pow(values[i] - center, order);
        weight_sum += static_cast<long double>(w);
    }
    if (weight_sum <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return static_cast<double>(normalize_weights ? weighted_sum / weight_sum : weighted_sum);
}

double weighted_univariate_cumulant(
    const double* values,
    const double* weights,
    int sample_count,
    int order,
    bool normalize_weights,
    bool zero_first_order_moment) {
    if (order <= 0) {
        return 0.0;
    }
    std::vector<double> moments(static_cast<size_t>(order) + 1, 0.0);
    std::vector<double> cumulants(static_cast<size_t>(order) + 1, 0.0);
    for (int n = 1; n <= order; ++n) {
        moments[static_cast<size_t>(n)] = weighted_univariate_moment(
            values, weights, sample_count, n, normalize_weights, zero_first_order_moment);
    }
    cumulants[1] = moments[1];
    for (int n = 2; n <= order; ++n) {
        long double value = moments[static_cast<size_t>(n)];
        for (int j = 1; j < n; ++j) {
            value -= static_cast<long double>(binomial(n - 1, j - 1)) *
                static_cast<long double>(cumulants[static_cast<size_t>(j)]) *
                static_cast<long double>(moments[static_cast<size_t>(n - j)]);
        }
        cumulants[static_cast<size_t>(n)] = static_cast<double>(value);
    }
    return cumulants[static_cast<size_t>(order)];
}

} // namespace weighted_moment_core
