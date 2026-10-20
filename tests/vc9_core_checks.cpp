#include "../shared/WeightedMomentCore.h"

#include <cmath>
#include <iostream>
#include <vector>

namespace wmc = weighted_moment_core;

static bool close_enough(double actual, double expected, double tolerance) {
    return std::fabs(actual - expected) <= tolerance;
}

static int fail(const char* label, double actual, double expected) {
    std::cerr << label << ": actual=" << actual << " expected=" << expected << std::endl;
    return 1;
}

static double central_moment(const double* values, int count, int order) {
    int index = 0;
    long double sum = 0.0;
    long double mean = 0.0;
    for(index = 0; index < count; ++index) mean += values[index];
    mean /= (long double)count;
    for(index = 0; index < count; ++index) {
        long double power = 1.0;
        int exponent = 0;
        const long double delta = (long double)values[index] - mean;
        for(exponent = 0; exponent < order; ++exponent) power *= delta;
        sum += power;
    }
    return (double)(sum / (long double)count);
}

int main() {
    const double tolerance = 1.0e-12;
    const double nonzero_mean[] = {0.0, 1.0, 3.0};
    const double coordinate_values[] = {0.1, 0.2, 0.3};
    const double coordinate_x[] = {-1.0, 0.0, 1.0};
    const double coordinate_y[] = {0.0, 0.0, 0.0};
    const double coordinate_z[] = {0.0, 0.0, 0.0};
    std::vector<double> weights;
    wmc::WeightDiagnostics diagnostics;
    wmc::WeightedMomentOptions options;

    if(!close_enough(wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 1, true, false), 4.0 / 3.0, tolerance)) return fail("kappa1", wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 1, true, false), 4.0 / 3.0);
    if(!close_enough(wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 2, true, false), 14.0 / 9.0, tolerance)) return fail("kappa2", wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 2, true, false), 14.0 / 9.0);
    if(!close_enough(wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 3, true, false), 20.0 / 27.0, tolerance)) return fail("kappa3", wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 3, true, false), 20.0 / 27.0);
    if(!close_enough(wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 4, true, false), -98.0 / 27.0, tolerance)) return fail("kappa4", wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 4, true, false), -98.0 / 27.0);

    {
        const double shifted_values[] = {1000000.0, 1000001.0, 1000003.0};
        const double m2 = central_moment(nonzero_mean, 3, 2);
        const double m3 = central_moment(nonzero_mean, 3, 3);
        const double m4 = central_moment(nonzero_mean, 3, 4);
        const double m5 = central_moment(nonzero_mean, 3, 5);
        const double m6 = central_moment(nonzero_mean, 3, 6);
        const double independent_k4 = m4 - 3.0 * m2 * m2;
        const double independent_k5 = m5 - 10.0 * m3 * m2;
        const double independent_k6 = m6 - 15.0 * m4 * m2 - 10.0 * m3 * m3 + 30.0 * m2 * m2 * m2;
        int order = 0;
        for(order = 2; order <= 8; ++order) {
            const double centered = wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, order, true, true);
            const double shifted = wmc::weighted_univariate_cumulant(shifted_values, NULL, 3, order, true, true);
            if(!close_enough(centered, shifted, 1.0e-8)) return fail("cumulant shift invariance", shifted, centered);
        }
        if(!close_enough(wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 4, true, true), independent_k4, tolerance)) return fail("independent kappa4", wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 4, true, true), independent_k4);
        if(!close_enough(wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 5, true, true), independent_k5, tolerance)) return fail("independent kappa5", wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 5, true, true), independent_k5);
        if(!close_enough(wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 6, true, true), independent_k6, tolerance)) return fail("independent kappa6", wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 6, true, true), independent_k6);
        if(!close_enough(independent_k4 / (m2 * m2), -1.5, tolerance)) return fail("lambda4 reference", independent_k4 / (m2 * m2), -1.5);
        if(!close_enough(
                wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 5, true, true) / std::pow(m2, 2.5),
                independent_k5 / std::pow(m2, 2.5), tolerance)) return fail("lambda5 reference", wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 5, true, true) / std::pow(m2, 2.5), independent_k5 / std::pow(m2, 2.5));
        if(!close_enough(
                wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 6, true, true) / std::pow(m2, 3.0),
                independent_k6 / std::pow(m2, 3.0), tolerance)) return fail("lambda6 reference", wmc::weighted_univariate_cumulant(nonzero_mean, NULL, 3, 6, true, true) / std::pow(m2, 3.0), independent_k6 / std::pow(m2, 3.0));
    }

    weights.assign(3, 1.0);
    diagnostics = wmc::compute_weight_diagnostics(&weights[0], 3);
    if(!close_enough(diagnostics.effective_sample_size, 3.0, tolerance)) return fail("uniform ESS", diagnostics.effective_sample_size, 3.0);

    options.mode = wmc::WeightMode_SoftHighPorosity;
    options.high_porosity_threshold = 0.25;
    options.soft_transition_width = 0.02;
    options.soft_weight_strength = 0.8;
    if(!close_enough(wmc::compute_weight(wmc::WeightedSample(), options), 0.0, tolerance)) return fail("soft missing sample", wmc::compute_weight(wmc::WeightedSample(), options), 0.0);
    weights = wmc::compute_univariate_weights(coordinate_values, 3, coordinate_x, coordinate_y, coordinate_z, options);
    if(!close_enough(weights[1], 0.260686544016995, 1.0e-12)) return fail("soft mid-range", weights[1], 0.260686544016995);
    options.soft_weight_strength = 0.0;
    weights = wmc::compute_univariate_weights(coordinate_values, 3, coordinate_x, coordinate_y, coordinate_z, options);
    if(!close_enough(weights[0], 1.0, tolerance) || !close_enough(weights[1], 1.0, tolerance) || !close_enough(weights[2], 1.0, tolerance)) return fail("soft eta zero", weights[1], 1.0);

    options = wmc::WeightedMomentOptions();
    options.mode = wmc::WeightMode_CoordinateGaussian;
    options.bandwidth = 1.0;
    options.coordinate_scale_x = 1.0;
    options.coordinate_scale_y = 0.0;
    options.coordinate_scale_z = 0.0;
    weights = wmc::compute_univariate_weights(coordinate_values, 3, coordinate_x, coordinate_y, coordinate_z, options);
    if(!close_enough(weights[1], 1.0, tolerance)) return fail("coordinate center", weights[1], 1.0);
    if(!close_enough(weights[0], std::exp(-0.5), tolerance)) return fail("coordinate one scale", weights[0], std::exp(-0.5));

    options = wmc::WeightedMomentOptions();
    options.mode = wmc::WeightMode_CoordinateGaussian;
    options.bandwidth = 1.0;
    if(!wmc::configure_relative_window_coordinate_gaussian(
        &options, 10.0, 20.0, 30.0, 2.0, 4.0, 8.0, 0, 0, 0, 3, 3, 1)) {
        std::cerr << "relative coordinate setup failed" << std::endl;
        return 1;
    }
    if(!close_enough(options.center_x, 12.0, tolerance) ||
       !close_enough(options.center_y, 24.0, tolerance) ||
       !close_enough(options.center_z, 30.0, tolerance) ||
       !close_enough(options.coordinate_scale_x, 2.0, tolerance) ||
       !close_enough(options.coordinate_scale_y, 4.0, tolerance) ||
       !close_enough(options.coordinate_scale_z, 0.0, tolerance)) {
        std::cerr << "relative coordinate geometry mismatch" << std::endl;
        return 1;
    }
    {
        wmc::WeightedSample left;
        wmc::WeightedSample center;
        wmc::WeightedSample right;
        left.x = 10.0; left.y = 24.0; left.z = 30.0;
        center.x = 12.0; center.y = 24.0; center.z = 30.0;
        right.x = 14.0; right.y = 24.0; right.z = 30.0;
        if(!close_enough(wmc::compute_weight(center, options), 1.0, tolerance)) return fail("relative coordinate center", wmc::compute_weight(center, options), 1.0);
        if(!close_enough(wmc::compute_weight(left, options), std::exp(-0.5), tolerance)) return fail("relative coordinate left", wmc::compute_weight(left, options), std::exp(-0.5));
        if(!close_enough(wmc::compute_weight(right, options), std::exp(-0.5), tolerance)) return fail("relative coordinate right", wmc::compute_weight(right, options), std::exp(-0.5));
    }

    std::cout << "VC9 shared weighted core checks passed: " << wmc::weighting_formula_version() << std::endl;
    return 0;
}
