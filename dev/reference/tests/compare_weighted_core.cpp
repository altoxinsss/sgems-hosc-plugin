#include "WeightedMomentCore.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace wmc = weighted_moment_core;

namespace {

bool close_enough(double actual, double expected, double tolerance = 1.0e-12) {
    return std::abs(actual - expected) <= tolerance * std::max(1.0, std::abs(expected));
}

int require(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAILED: " << message << "\n";
        return 1;
    }
    return 0;
}

template <typename Action>
bool rejects_invalid_weight(Action action) {
    try {
        action();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

}

int main() {
    int failures = 0;
    const double values[] = {1.0, 2.0, 4.0, 7.0};
    const double weights[] = {1.0, 2.0, 3.0, 4.0};

    const double centered_moment4 = wmc::weighted_univariate_moment(
        values, weights, 4, 4, true, true);
    const double centered_cumulant4 = wmc::weighted_univariate_cumulant(
        values, weights, 4, 4, true, true);
    const double raw_cumulant4 = wmc::weighted_univariate_cumulant(
        values, weights, 4, 4, true, false);
    const double centered_cumulant2 = wmc::weighted_univariate_cumulant(
        values, weights, 4, 2, true, true);
    const double centered_cumulant3 = wmc::weighted_univariate_cumulant(
        values, weights, 4, 3, true, true);
    const double centered_cumulant5 = wmc::weighted_univariate_cumulant(
        values, weights, 4, 5, true, true);
    const double centered_cumulant6 = wmc::weighted_univariate_cumulant(
        values, weights, 4, 6, true, true);
    // Independently evaluated weighted central moments for values={1,2,4,7}
    // and weights={1,2,3,4}; mu=4.5, m2=5.05, m3=-1.2, m4=38.4625.
    const double expected_kappa2 = 5.05;
    const double expected_kappa3 = -1.2;
    const double expected_kappa4 = 38.4625 - 3.0 * expected_kappa2 * expected_kappa2;
    // The same independently evaluated weighted central moments give
    // m5=-33 and m6=330.315625.  These identities exercise the recurrence
    // used by the formal fifth- and sixth-order feature maps.
    const double expected_kappa5 = -33.0 - 10.0 * expected_kappa3 * expected_kappa2;
    const double expected_kappa6 = 330.315625 - 15.0 * 38.4625 * expected_kappa2
        - 10.0 * expected_kappa3 * expected_kappa3 + 30.0 * expected_kappa2 * expected_kappa2 * expected_kappa2;
    failures += require(close_enough(centered_cumulant2, expected_kappa2),
        "weighted centered kappa2 must equal the explicit second central moment");
    failures += require(close_enough(centered_cumulant3, expected_kappa3),
        "weighted centered kappa3 must equal the explicit third central moment");
    failures += require(close_enough(centered_cumulant4, expected_kappa4),
        "weighted centered kappa4 must equal m4 minus three times m2 squared");
    failures += require(close_enough(centered_cumulant5, expected_kappa5),
        "weighted centered kappa5 must match the explicit fifth-order cumulant identity");
    failures += require(close_enough(centered_cumulant6, expected_kappa6),
        "weighted centered kappa6 must match the explicit sixth-order cumulant identity");
    failures += require(close_enough(raw_cumulant4, expected_kappa4),
        "raw nonzero-mean recurrence must recover the fourth cumulant");
    failures += require(close_enough(wmc::weighted_univariate_cumulant(values, weights, 4, 2), expected_kappa2),
        "raw nonzero-mean recurrence must recover the variance");
    failures += require(close_enough(wmc::weighted_univariate_cumulant(values, weights, 4, 3), expected_kappa3),
        "raw nonzero-mean recurrence must recover the third cumulant");

    const double bivariate_values[][2] = {{1.0, 2.0}, {2.0, 4.0}, {4.0, 7.0}, {7.0, 11.0}};
    wmc::WeightedSample bivariate_samples[4];
    for (int i = 0; i < 4; ++i) {
        bivariate_samples[i].values = bivariate_values[i];
        bivariate_samples[i].value_count = 2;
    }
    const int bivariate_multiplicities[] = {2, 2};
    const wmc::PreparedCumulantPlan prepared = wmc::prepare_cumulant_plan(bivariate_multiplicities, 2);
    const double planned_bivariate = wmc::evaluate_weighted_cumulant(
        bivariate_samples, 4, prepared, wmc::WeightedMomentOptions{});
    const double legacy_bivariate = wmc::weighted_cumulant_from_multiplicities(
        bivariate_samples, 4, bivariate_multiplicities, 2, wmc::WeightedMomentOptions{});
    failures += require(prepared.order == 4 && !prepared.terms.empty() && !prepared.required_masks.empty(),
        "prepared cumulant plan must retain fourth-order partition structure");
    failures += require(close_enough(planned_bivariate, legacy_bivariate),
        "prepared and compatibility cumulant evaluation must agree");
    const wmc::WeightDiagnostics diagnostics = wmc::compute_weight_diagnostics(
        bivariate_samples, 4, wmc::WeightedMomentOptions{});
    failures += require(diagnostics.positive_weight_count == 4 && close_enough(diagnostics.effective_sample_size, 4.0),
        "uniform weight diagnostics must report the full effective sample size");

    wmc::WeightedMomentOptions external;
    external.mode = wmc::WeightMode::External;
    wmc::WeightedSample scalar_samples[4];
    double scaled_weights[4];
    for (int i = 0; i < 4; ++i) {
        scalar_samples[i].values = values + i;
        scalar_samples[i].value_count = 1;
        scalar_samples[i].weight = weights[i];
        scaled_weights[i] = weights[i] * 17.0;
    }
    const int fourth[] = {4};
    const int first[] = {1};
    failures += require(close_enough(wmc::weighted_moment(scalar_samples, 4, first, 1, external), 4.5),
        "external weights must give the explicit weighted mean");
    failures += require(close_enough(wmc::weighted_cumulant_from_multiplicities(
        scalar_samples, 4, fourth, 1, external), expected_kappa4),
        "external scalar Bell cumulant must equal the closed-form result");
    failures += require(close_enough(wmc::weighted_univariate_cumulant(
        values, scaled_weights, 4, 4), expected_kappa4),
        "normalized cumulants must be invariant to positive weight scaling");
    const auto external_diagnostics = wmc::compute_weight_diagnostics(scalar_samples, 4, external);
    failures += require(close_enough(external_diagnostics.weight_sum, 10.0) &&
        close_enough(external_diagnostics.effective_sample_size, 100.0 / 30.0),
        "external weight diagnostics must match their closed forms");
    for (int i = 0; i < 4; ++i) {
        scalar_samples[i].weight = scaled_weights[i];
    }
    failures += require(close_enough(wmc::weighted_cumulant_from_multiplicities(
        scalar_samples, 4, fourth, 1, external), expected_kappa4),
        "joint API normalization must also be invariant to weight scaling");

    // Repeated variables use the complete-tuple weight once, not one weight
    // factor per occurrence. For this weighted sample, central mu22=83.0125,
    // mu20=5.05, mu02=10.45 and mu11=7.25.
    for (int i = 0; i < 4; ++i) {
        bivariate_samples[i].weight = weights[i];
    }
    const double joint_expected = 83.0125 - 5.05 * 10.45 - 2.0 * 7.25 * 7.25;
    failures += require(close_enough(wmc::evaluate_weighted_cumulant(
        bivariate_samples, 4, prepared, external), joint_expected),
        "joint external kappa22 must match the explicit central-moment identity");
    external.zero_first_order_moment = true;
    failures += require(close_enough(wmc::evaluate_weighted_cumulant(
        bivariate_samples, 4, prepared, external), joint_expected),
        "joint external kappa22 must agree after centering");

    for (int i = 0; i < 4; ++i) {
        bivariate_samples[i].weight = 7.0;
    }
    failures += require(close_enough(wmc::evaluate_weighted_cumulant(
        bivariate_samples, 4, prepared, external), planned_bivariate),
        "equal external weights must reproduce uniform joint cumulants");
    const double invalid_weights[] = {-1.0, std::numeric_limits<double>::infinity(),
        std::numeric_limits<double>::quiet_NaN()};
    for (double invalid : invalid_weights) {
        scalar_samples[0].weight = invalid;
        double bad_array[] = {invalid, 2.0, 3.0, 4.0};
        failures += require(rejects_invalid_weight([&] {
            wmc::weighted_moment(scalar_samples, 4, first, 1, external);
        }), "generic sample API must reject negative or non-finite weights");
        failures += require(rejects_invalid_weight([&] {
            wmc::weighted_univariate_cumulant(values, bad_array, 4, 4);
        }), "univariate array API must reject negative or non-finite weights");
    }
    const double zero_weights[] = {0.0, 0.0, 0.0, 0.0};
    failures += require(std::isnan(wmc::weighted_univariate_cumulant(values, zero_weights, 4, 4)),
        "an empty positive-weight support must produce an undefined statistic");
    for (int i = 0; i < 4; ++i) {
        scalar_samples[i].weight = 0.0;
    }
    failures += require(std::isnan(wmc::weighted_cumulant_from_multiplicities(
        scalar_samples, 4, fourth, 1, external)),
        "joint API must also report undefined statistics for zero total weight");
    const double ignored_value[] = {std::numeric_limits<double>::quiet_NaN(), 2.0};
    const double ignored_weights[] = {0.0, 1.0};
    failures += require(close_enough(wmc::weighted_univariate_cumulant(
        ignored_value, ignored_weights, 2, 2), 0.0),
        "a zero-weight observation must not contaminate the statistic");
    failures += require(close_enough(wmc::weighted_univariate_moment(values, weights, 4, 1, false), 45.0),
        "the generic raw weighted-sum moment API must remain available");
    wmc::WeightedMomentOptions unnormalized = external;
    unnormalized.normalize_weights = false;
    failures += require(rejects_invalid_weight([&] {
        wmc::evaluate_weighted_cumulant(bivariate_samples, 4, prepared, unnormalized);
    }), "joint cumulants must not silently use unnormalized sums as moments");
    failures += require(rejects_invalid_weight([&] {
        wmc::weighted_univariate_cumulant(values, weights, 4, 4, false);
    }), "univariate cumulants must not silently use unnormalized sums as moments");
    failures += require(rejects_invalid_weight([&] {
        wmc::compute_univariate_weights(values, 4, external);
    }), "an external mode without supplied weights must not fall back to uniform");

    std::cout << std::setprecision(17)
              << "centered_moment4=" << centered_moment4 << "\n"
              << "centered_cumulant2=" << centered_cumulant2 << "\n"
              << "centered_cumulant3=" << centered_cumulant3 << "\n"
              << "centered_cumulant4=" << centered_cumulant4 << "\n"
              << "centered_cumulant5=" << centered_cumulant5 << "\n"
              << "centered_cumulant6=" << centered_cumulant6 << "\n"
              << "raw_cumulant4=" << raw_cumulant4 << "\n"
              << "weighting_formula_version=" << wmc::weighting_formula_version() << "\n"
              << "checks=" << (failures == 0 ? "PASS" : "FAIL") << "\n";
    return failures == 0 ? 0 : 1;
}
