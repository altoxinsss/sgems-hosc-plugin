#ifndef WEIGHTED_MOMENT_CORE_H
#define WEIGHTED_MOMENT_CORE_H

#include "NativeTypes.h"
#include <cstddef>
#include <string>
#include <vector>

namespace weighted_moment_core {

enum WeightMode {
    WeightMode_Uniform,
    WeightMode_External,
};

struct WeightedMomentOptions {
    WeightMode mode;
    // False requests raw weighted sums from moment APIs. Cumulant APIs and
    // statistical CLI calculations require normalized moments.
    bool normalize_weights;
    bool zero_first_order_moment;
    WeightedMomentOptions():mode(WeightMode_Uniform),normalize_weights(true),zero_first_order_moment(false) {}
};

struct WeightedSample {
    const double* values;
    int value_count;
    double x;
    double y;
    double z;
    // A caller-supplied sample/complete-tuple weight. Multiplicities never
    // change this weight. Ignored in Uniform mode.
    double weight;
    WeightedSample():values(NULL),value_count(0),x(0.0),y(0.0),z(0.0),weight(1.0) {}
};

struct CumulantPartitionTerm {
    long double coefficient;
    std::vector<unsigned int> block_masks;
    CumulantPartitionTerm():coefficient(0.0) {}
};

// The combinatorial structure of a joint cumulant depends only on the
// multiplicity vector.  It can therefore be compiled once and reused for
// every local-template replicate set.
struct PreparedCumulantPlan {
    std::vector<int> labels;
    std::vector<unsigned int> required_masks;
    std::vector<CumulantPartitionTerm> terms;
    int order;
    PreparedCumulantPlan():order(0) {}
};

struct WeightDiagnostics {
    int positive_weight_count;
    double weight_sum;
    double weight_square_sum;
    double effective_sample_size;
    WeightDiagnostics():positive_weight_count(0),weight_sum(0.0),weight_square_sum(0.0),effective_sample_size(0.0) {}
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
    WeightedMomentOptions* resolved_options = NULL);

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
