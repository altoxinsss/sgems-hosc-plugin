#include "WeightedMomentCore.h"

#include <cmath>
#include <iomanip>
#include <iostream>
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

    const double tuple_values[] = {0.26, 0.34};
    wmc::WeightedSample tuple;
    tuple.values = tuple_values;
    tuple.value_count = 2;
    wmc::WeightedMomentOptions fixed;
    fixed.mode = wmc::WeightMode::ValueGaussian;
    fixed.reference_value = 0.26;
    fixed.bandwidth = 0.08;
    failures += require(close_enough(wmc::compute_weight(tuple, fixed), 1.0),
        "legacy fixed Gaussian must use the first component");
    fixed.joint_geometric_mean = true;
    failures += require(close_enough(wmc::compute_weight(tuple, fixed), std::exp(-0.25)),
        "joint fixed Gaussian must use the geometric mean of component weights");

    wmc::WeightedMomentOptions soft;
    soft.mode = wmc::WeightMode::SoftHighPorosity;
    soft.high_porosity_threshold = 0.25;
    soft.soft_transition_width = 0.02;
    soft.soft_weight_strength = 1.0;
    const double soft_values[] = {0.20, 0.25, 0.30};
    const std::vector<double> soft_weights = wmc::compute_univariate_weights(soft_values, 3, soft);
    failures += require(close_enough(soft_weights[1], 0.5),
        "soft high-porosity weight must equal 0.5 at the threshold");
    failures += require(soft_weights[0] < soft_weights[1] && soft_weights[1] < soft_weights[2],
        "soft high-porosity weights must increase with value");

    wmc::WeightedMomentOptions coordinate;
    coordinate.mode = wmc::WeightMode::CoordinateGaussian;
    coordinate.center_x = 10.0;
    coordinate.center_y = 20.0;
    coordinate.center_z = 30.0;
    coordinate.bandwidth = 2.0;
    wmc::WeightedSample coordinate_center;
    coordinate_center.x = 10.0;
    coordinate_center.y = 20.0;
    coordinate_center.z = 30.0;
    wmc::WeightedSample coordinate_one_bandwidth;
    coordinate_one_bandwidth.x = 12.0;
    coordinate_one_bandwidth.y = 20.0;
    coordinate_one_bandwidth.z = 30.0;
    failures += require(close_enough(wmc::compute_weight(coordinate_center, coordinate), 1.0),
        "coordinate Gaussian must equal one at its centre");
    failures += require(close_enough(wmc::compute_weight(coordinate_one_bandwidth, coordinate), std::exp(-0.5)),
        "coordinate Gaussian must equal exp(-1/2) one bandwidth from its centre");

    // Formal moving-window experiments use relative coordinates: an edge of a
    // 3-cell axis is one half-width from the centre, irrespective of spacing.
    coordinate.coordinate_scale_x = 2.0;
    coordinate.coordinate_scale_y = 3.0;
    coordinate.coordinate_scale_z = 4.0;
    coordinate_one_bandwidth.x = 12.0;
    failures += require(close_enough(wmc::compute_weight(coordinate_one_bandwidth, coordinate), std::exp(-0.125)),
        "relative-coordinate Gaussian must use the configured axis half-width");

    wmc::WeightedMomentOptions adaptive;
    adaptive.mode = wmc::WeightMode::AdaptiveUpperTail;
    adaptive.adaptive_global_threshold = 0.25;
    adaptive.adaptive_mean_offset = 0.0;
    adaptive.adaptive_transition_scale = 1.0;
    adaptive.adaptive_std_floor = 1.0e-6;
    wmc::WeightedMomentOptions resolved;
    const std::vector<double> adaptive_weights = wmc::compute_univariate_weights(
        soft_values, 3, adaptive, &resolved);
    failures += require(resolved.adaptive_context_ready, "adaptive context must be recorded");
    failures += require(close_enough(resolved.adaptive_context_threshold, 0.25),
        "adaptive threshold must respect the global floor");
    failures += require(adaptive_weights[0] < adaptive_weights[1] && adaptive_weights[1] < adaptive_weights[2],
        "adaptive upper-tail weights must increase with value");
    failures += require(wmc::weight_scheme_id(fixed) ==
        "fixed_value_gaussian.v1+joint_geometric_mean.v1",
        "scheme identifier must encode joint aggregation");
    failures += require(std::string(wmc::weighting_formula_version()) == "2026-10-10.v2",
        "relative-window coordinate Gaussian must report formula version 2026-10-10.v2");

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
