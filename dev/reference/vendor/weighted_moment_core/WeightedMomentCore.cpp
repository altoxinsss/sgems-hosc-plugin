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

const char* const kWeightingFormulaVersion = "2026-11-06.generic.v1";

double validated_weight(double weight) {
    if (!std::isfinite(weight) || weight < 0.0) {
        throw std::invalid_argument("External weights must be finite and non-negative.");
    }
    return weight;
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
    switch (options.mode) {
    case WeightMode::Uniform:
        return 1.0;
    case WeightMode::External:
        return validated_weight(sample.weight);
    default:
        throw std::invalid_argument("Unknown weight mode.");
    }
}

WeightedMomentOptions prepare_weight_context(
    const WeightedSample* samples,
    int sample_count,
    const WeightedMomentOptions& options) {
    // Compatibility hook: supplied weights do not depend on the values or
    // coordinates of a support and require no fitted context.
    (void)samples;
    (void)sample_count;
    return options;
}

std::vector<double> compute_univariate_weights(
    const double* values,
    int sample_count,
    const WeightedMomentOptions& options,
    WeightedMomentOptions* resolved_options) {
    if (values == nullptr || sample_count <= 0) {
        throw std::invalid_argument("compute_univariate_weights received invalid input.");
    }
    if (options.mode != WeightMode::Uniform) {
        throw std::invalid_argument(
            "compute_univariate_weights only generates uniform weights; supply external weights explicitly.");
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
    case WeightMode::External:
        return "external";
    case WeightMode::Uniform:
        return "uniform";
    default:
        throw std::invalid_argument("Unknown weight mode.");
    }
}

std::string weight_scheme_id(const WeightedMomentOptions& options) {
    return std::string(weight_mode_name(options.mode)) + ".v1";
}

std::string weight_options_manifest(const WeightedMomentOptions& options) {
    std::ostringstream stream;
    stream.precision(17);
    stream << "weighting_formula_version=" << weighting_formula_version() << "\n";
    stream << "weighting_scheme_id=" << weight_scheme_id(options) << "\n";
    stream << "weight_mode=" << weight_mode_name(options.mode) << "\n";
    stream << "normalize_weights=" << (options.normalize_weights ? "true" : "false") << "\n";
    stream << "zero_first_order_moment=" << (options.zero_first_order_moment ? "true" : "false") << "\n";
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
    if (!options.normalize_weights) {
        throw std::invalid_argument("Cumulants require normalized moments; normalize_weights must be true.");
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
    if (labels.size() >= sizeof(int) * 8U - 1U) {
        throw std::runtime_error("Joint-moment order exceeds the existing partition-mask capacity.");
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
    if (weights != nullptr) {
        for (int i = 0; i < sample_count; ++i) {
            validated_weight(weights[i]);
        }
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
    if (!normalize_weights) {
        throw std::invalid_argument("Cumulants require normalized moments; normalize_weights must be true.");
    }
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
