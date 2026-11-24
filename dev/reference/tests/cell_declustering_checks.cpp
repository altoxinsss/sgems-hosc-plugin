#include "cell_declustering.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>

using namespace without_sgems;

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

void close(double actual, double expected, const std::string& message) {
    require(std::abs(actual - expected) <= 2e-14 * std::max(1.0, std::abs(expected)), message);
}

template<class Operation> void rejects(Operation operation, const std::string& message) {
    try { operation(); } catch (const std::exception&) { return; }
    throw std::runtime_error(message);
}

PointSample point(int id, double x, double y, double z, double value, bool informed = true) {
    return {id, x, y, z, value, informed};
}

void check_weights(const CellDeclusteringResult& result, const std::vector<double>& expected,
    const std::string& name) {
    require(result.weights.size() == expected.size(), name + ": row count differs");
    for (std::size_t row = 0; row < expected.size(); ++row) {
        close(result.weights[row], expected[row], name + ": incorrect weight at row " + std::to_string(row));
    }
    close(std::accumulate(result.weights.begin(), result.weights.end(), 0.0), 1.0, name + ": weights do not sum to one");
}

}

int main() {
    try {
        CellDeclusteringOptions fixed;
        fixed.origin_at_data_min = false;
        fixed.shifts_per_axis = 1;

        // Two occupied cells, with occupancies 2 and 1. Missing rows are not
        // cell occupants even if their coordinates are outside the valid cloud.
        PointSetData negative;
        negative.dimensions = 2;
        negative.nodes = {
            point(1, -1.875, .25, 0, 1), point(2, -1.125, .25, 0, 3),
            point(3, .25, .25, 0, 8), point(4, -1000, -1000, 0, 100, false)
        };
        auto result = calculate_cell_declustering(negative, fixed);
        check_weights(result, {.25, .25, .5, 0}, "negative and missing fixture");
        close(result.weighted_mean, 5, "Weighted mean incorrect");
        require(result.valid_count == 3 && result.origin_count == 1 &&
            result.min_occupied_cells == 2 && result.max_occupied_cells == 2, "Occupancy diagnostics incorrect");

        // The two x origins produce (1/4,1/4,1/2) and (1/2,1/4,1/4).
        // Both y origins produce the same result, so their average is exact.
        PointSetData shifted;
        shifted.dimensions = 2;
        shifted.nodes = {point(1, .125, .25, 0, 2), point(2, .625, .25, 0, 4), point(3, 1.125, .25, 0, 8)};
        auto offsets = fixed;
        offsets.shifts_per_axis = 2;
        result = calculate_cell_declustering(shifted, offsets);
        check_weights(result, {.375, .25, .375}, "fresh occupancy across offsets");
        require(result.origin_count == 4, "Two-dimensional origin count incorrect");
        close(result.weighted_mean, 4.75, "Offset-averaged mean incorrect");
        require(calculate_cell_declustering(shifted, offsets).weights == result.weights,
            "Repeated calculation is not deterministic");
        auto variable_occupancy = shifted;
        variable_occupancy.dimensions = 1;
        variable_occupancy.nodes[2].x = 1.625;
        result = calculate_cell_declustering(variable_occupancy, offsets);
        // Origin 0 has two occupied cells; origin 1/2 has three. Normalize
        // within each origin before averaging, rather than pooling counts.
        check_weights(result, {7.0/24,7.0/24,5.0/12}, "separately normalized origin lattices");
        require(result.min_occupied_cells == 2 && result.max_occupied_cells == 3,
            "Varying occupancy diagnostics incorrect");
        variable_occupancy.nodes[0].value = -1700;
        variable_occupancy.nodes[1].value = 450000;
        check_weights(calculate_cell_declustering(variable_occupancy, offsets), result.weights,
            "Weights depend on property values rather than sampling geometry");

        // Each z cell must use the z width, independently of the x width.
        PointSetData layers;
        layers.dimensions = 3;
        layers.nodes = {point(1, .125, .125, .125, 0), point(2, .375, .125, .375, 2), point(3, .625, .125, 1.125, 8)};
        auto anisotropic = fixed;
        anisotropic.cell_size = {{10, 10, 1}};
        check_weights(calculate_cell_declustering(layers, anisotropic), {.25, .25, .5}, "independent z cell width");
        anisotropic.cell_size[2] = 10;
        check_weights(calculate_cell_declustering(layers, anisotropic), {1.0/3, 1.0/3, 1.0/3}, "wide z cell");
        anisotropic.cell_size[0] = .25;
        check_weights(calculate_cell_declustering(layers, anisotropic), {1.0/3, 1.0/3, 1.0/3}, "independent x cell width");

        // Coincident observations retain separate rows, sharing their cell mass.
        PointSetData duplicate;
        duplicate.dimensions = 1;
        duplicate.nodes = {point(1, .125, 0, 0, 1), point(2, .125, 0, 0, 9), point(3, 2.125, 0, 0, 5)};
        check_weights(calculate_cell_declustering(duplicate, fixed), {.25, .25, .5}, "duplicate observations");
        auto default_options = CellDeclusteringOptions{};
        // XY declustering of equally sampled vertical wells must retain equal
        // within-well layer weights despite unequal Z coordinates.
        PointSetData wells;
        wells.dimensions=3;
        wells.nodes={point(1,.125,.25,1,2),point(2,.125,.25,2,3),
            point(3,.625,.25,1,4),point(4,.625,.25,2,5),
            point(5,1.125,.25,1,8),point(6,1.125,.25,2,9)};
        auto horizontal=fixed;
        horizontal.dimensions=2;
        horizontal.origin_offsets={{{0,0,0}},{{.5,0,0}}};
        check_weights(calculate_cell_declustering(wells,horizontal),
            {.1875,.1875,.125,.125,.1875,.1875},"explicit XY origins with repeated well layers");
        horizontal.origin_offsets[1][2]=.5;
        rejects([&]{calculate_cell_declustering(wells,horizontal);},"inactive fractional origin accepted");
        horizontal.origin_offsets[1][2]=0;
        horizontal.origin_offsets[1][0]=1;
        rejects([&]{calculate_cell_declustering(wells,horizontal);},"out-of-cell fractional origin accepted");
        const auto data_min = calculate_cell_declustering(negative, default_options);
        require(data_min.origin == std::array<double,3>{{-1.875,.25,0}}, "Missing row affected data-min origin");

        // Translation, common scale changes, and row permutation preserve the
        // physical cell membership when the origin convention is transformed too.
        auto reference = calculate_cell_declustering(layers, default_options);
        auto translated = layers;
        for (auto& sample : translated.nodes) { sample.x += 16; sample.y -= 32; sample.z += 8; }
        check_weights(calculate_cell_declustering(translated, default_options), reference.weights, "data-min translation invariance");
        auto scaled = layers;
        for (auto& sample : scaled.nodes) { sample.x *= 8; sample.y *= 8; sample.z *= 8; }
        auto scaled_options = default_options;
        for (auto& width : scaled_options.cell_size) width *= 8;
        check_weights(calculate_cell_declustering(scaled, scaled_options), reference.weights, "physical scale invariance");
        std::reverse(scaled.nodes.begin(), scaled.nodes.end());
        auto reverse_expected = reference.weights;
        std::reverse(reverse_expected.begin(), reverse_expected.end());
        check_weights(calculate_cell_declustering(scaled, scaled_options), reverse_expected, "row permutation invariance");
        auto explicit_translated = fixed;
        explicit_translated.origin = {{16,-32,8}};
        check_weights(calculate_cell_declustering(translated, explicit_translated),
            calculate_cell_declustering(layers, fixed).weights, "explicit-origin translation invariance");

        // Boundary cells are half-open; floor, not truncation toward zero,
        // determines negative-coordinate membership.
        PointSetData boundary;
        boundary.dimensions = 1;
        boundary.nodes = {point(1,-1,0,0,1),point(2,-.5,0,0,2),point(3,0,0,0,3),point(4,1,0,0,4)};
        check_weights(calculate_cell_declustering(boundary,fixed), {1.0/6,1.0/6,1.0/3,1.0/3}, "cell boundaries");

        auto one = layers;
        one.nodes.resize(1);
        result = calculate_cell_declustering(one, default_options);
        check_weights(result, {1}, "single valid sample");
        require(result.origin_count == 27 && result.min_occupied_cells == 1 && result.max_occupied_cells == 1,
            "Three-dimensional origin diagnostics incorrect");
        auto maximum_origins = default_options;
        maximum_origins.shifts_per_axis = 16;
        require(calculate_cell_declustering(one, maximum_origins).origin_count == 4096,
            "Allowed origin-count boundary rejected");
        one.nodes[0].value = std::numeric_limits<double>::max();
        require(calculate_cell_declustering(one, default_options).weighted_mean == one.nodes[0].value,
            "Finite extreme value caused mean overflow");
        PointSetData extreme;
        extreme.dimensions = 1;
        const double largest = std::numeric_limits<double>::max();
        extreme.nodes = {point(1,-largest,0,0,-largest), point(2,largest,0,0,largest)};
        auto extreme_options = default_options;
        extreme_options.cell_size[0] = largest;
        const auto extreme_result = calculate_cell_declustering(extreme, extreme_options);
        check_weights(extreme_result, {.5,.5}, "extreme finite coordinates with large cells");
        close(extreme_result.weighted_mean, 0, "Extreme opposite values overflowed mean");

        for (double bad : {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
            auto invalid = fixed;
            invalid.cell_size[2] = bad;
            rejects([&] { calculate_cell_declustering(negative, invalid); }, "Invalid inactive-axis width accepted");
        }
        auto invalid = default_options;
        invalid.origin[1] = std::numeric_limits<double>::infinity();
        rejects([&] { calculate_cell_declustering(negative, invalid); }, "Nonfinite origin accepted");
        for (int bad : {0, -1, 17, std::numeric_limits<int>::max()}) {
            invalid = default_options;
            invalid.shifts_per_axis = bad;
            rejects([&] { calculate_cell_declustering(layers, invalid); }, "Invalid/excessive origin count accepted");
        }
        invalid = fixed;
        invalid.cell_size[0] = std::numeric_limits<double>::denorm_min();
        rejects([&] { calculate_cell_declustering(layers, invalid); }, "Overflowing cell index accepted");
        auto invalid_data = layers;
        invalid_data.nodes[0].value = std::numeric_limits<double>::quiet_NaN();
        rejects([&] { calculate_cell_declustering(invalid_data, fixed); }, "Informed NaN value accepted");
        invalid_data = layers;
        invalid_data.nodes[0].x = std::numeric_limits<double>::infinity();
        rejects([&] { calculate_cell_declustering(invalid_data, fixed); }, "Nonfinite informed coordinate accepted");
        invalid_data = layers;
        invalid_data.dimensions = 4;
        rejects([&] { calculate_cell_declustering(invalid_data, fixed); }, "Invalid coordinate dimension accepted");
        invalid_data = layers;
        for (auto& sample : invalid_data.nodes) sample.informed = false;
        rejects([&] { calculate_cell_declustering(invalid_data, fixed); }, "All-missing input accepted");
        invalid_data.nodes.clear();
        rejects([&] { calculate_cell_declustering(invalid_data, fixed); }, "Empty input accepted");

        std::cout << "Cell-declustering occupancy, offsets, physical coordinates, normalization, invariance, and validation checks passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
