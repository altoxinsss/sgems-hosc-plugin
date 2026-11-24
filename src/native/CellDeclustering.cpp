#include "NativeTypes.h"
#include <algorithm>
#include <cmath>

#include <limits>
#include <stdexcept>
#include <map>

namespace without_sgems {
namespace {

typedef Triple<__int64> Cell;

Vec3 coordinates(const PointSample& sample) {
    return Vec3(sample.x,sample.y,sample.z);
}

__int64 cell_index(double coordinate, double origin, double width,
    long double fractional_shift) {
    // Subtract the fractional origin shift after division, avoiding overflow
    // from forming origin + width * shift / shift_count in physical units.
    const long double delta = static_cast<long double>(coordinate) - static_cast<long double>(origin);
    // MSVC represents long double as double. Opposite extreme coordinates can
    // overflow their difference even when division by a large width is safe.
    const long double scaled_delta = without_sgems::finite(delta)
        ? delta / static_cast<long double>(width)
        : static_cast<long double>(coordinate) / static_cast<long double>(width)
            - static_cast<long double>(origin) / static_cast<long double>(width);
    const long double position = scaled_delta - fractional_shift;
    const long double index = std::floor(position);
    // Use an exclusive power-of-two upper bound. Converting INT64_MAX to
    // floating point can round it up to 2^63 on platforms with 64-bit doubles.
    const long double limit = std::ldexp(1.0L, 63);
    if (!without_sgems::finite(index) || index < -limit || index >= limit) {
        throw std::runtime_error("Cell declustering cell index exceeds the signed 64-bit range; enlarge cell_size or adjust origin.");
    }
    return static_cast<__int64>(index);
}

}

CellDeclusteringResult calculate_cell_declustering(const PointSetData& data,
    const CellDeclusteringOptions& options) {
    if (data.dimensions < 1 || data.dimensions > 3) {
        throw std::runtime_error("Cell declustering requires one, two, or three coordinate dimensions.");
    }
    const int dimensions = options.dimensions == 0 ? data.dimensions : options.dimensions;
    if (dimensions < 1 || dimensions > data.dimensions)
        throw std::runtime_error("Declustering dimensions must be between one and the input coordinate dimension.");
    for (int axis = 0; axis < 3; ++axis) {
        if (!without_sgems::finite(options.cell_size[axis]) || options.cell_size[axis] <= 0.0) {
            throw std::runtime_error("Cell declustering cell_size components must be finite and positive.");
        }
        if (!without_sgems::finite(options.origin[axis])) {
            throw std::runtime_error("Cell declustering origin components must be finite.");
        }
    }
    if (options.shifts_per_axis <= 0) {
        throw std::runtime_error("Cell declustering shifts_per_axis must be a positive integer.");
    }

    CellDeclusteringResult result;
    result.origin = options.origin;
    result.origin_count = 1;
    const std::size_t maximum_origins = 4096;
    for (int axis = 0; options.origin_offsets.empty() && axis < dimensions; ++axis) {
        const size_t count = static_cast<std::size_t>(options.shifts_per_axis);
        if (count > maximum_origins / result.origin_count) {
            throw std::runtime_error("Cell declustering Cartesian origin shifts exceed the limit of 4096; reduce shifts_per_axis.");
        }
        result.origin_count *= count;
    }
    if (!options.origin_offsets.empty()) {
        result.origin_count = options.origin_offsets.size();
        if (result.origin_count > maximum_origins)
            throw std::runtime_error("Explicit cell origins exceed the limit of 4096.");
        for (std::vector<Vec3>::const_iterator it_offset=(options.origin_offsets).begin();it_offset!=(options.origin_offsets).end();++it_offset) for (bool once_offset=true;once_offset;) for (const std::vector<Vec3>::value_type& offset=*it_offset;once_offset;once_offset=false) for (int axis = 0; axis < 3; ++axis)
            if (!without_sgems::finite(offset[axis]) || offset[axis] < 0 || offset[axis] >= 1 ||
                (axis >= dimensions && offset[axis] != 0))
                throw std::runtime_error("Origin offsets must be finite fractions in [0,1), with zero inactive components.");
    }

    std::vector<std::size_t> informed;
    informed.reserve(data.nodes.size());
    double largest_absolute_value = 0.0;
    for (std::size_t row = 0; row < data.nodes.size(); ++row) {
        const PointSample& sample = data.nodes[row];
        if (!sample.informed) continue;
        const Vec3 position = coordinates(sample);
        if (!without_sgems::finite(sample.value)) {
            throw std::runtime_error("Cell declustering informed sample values must be finite.");
        }
        for (int axis = 0; axis < 3; ++axis) {
            if (!without_sgems::finite(position[axis])) {
                throw std::runtime_error("Cell declustering informed sample coordinates must be finite.");
            }
            if (options.origin_at_data_min &&
                (informed.empty() || position[axis] < result.origin[axis])) {
                result.origin[axis] = position[axis];
            }
        }
        largest_absolute_value = std::max(largest_absolute_value, std::abs(sample.value));
        informed.push_back(row);
    }
    if (informed.empty()) {
        throw std::runtime_error("Cell declustering requires at least one informed sample.");
    }
    result.valid_count = informed.size();
    result.weights.assign(data.nodes.size(), 0.0);
    result.min_occupied_cells = informed.size();
    std::vector<long double> accumulated(informed.size(), 0.0L);
    std::vector<Cell> sample_cells(informed.size());

    for (std::size_t origin_number = 0; origin_number < result.origin_count; ++origin_number) {
        Triple<long double> shifts;
        std::size_t remaining = origin_number;
        for (int axis = 0; axis < dimensions; ++axis) {
            shifts[axis] = options.origin_offsets.empty()
                ? static_cast<long double>(remaining % static_cast<std::size_t>(options.shifts_per_axis)) / options.shifts_per_axis
                : options.origin_offsets[origin_number][axis];
            remaining /= static_cast<std::size_t>(options.shifts_per_axis);
        }
        // Occupancy must start empty for every origin.
        std::map<Cell, std::size_t> occupancy;
        for (std::size_t selected = 0; selected < informed.size(); ++selected) {
            const Vec3 position = coordinates(data.nodes[informed[selected]]);
            Cell cell;
            for (int axis = 0; axis < dimensions; ++axis) {
                cell[axis] = cell_index(position[axis], result.origin[axis], options.cell_size[axis],
                    shifts[axis]);
            }
            sample_cells[selected] = cell;
            ++occupancy[cell];
        }
        result.min_occupied_cells = std::min(result.min_occupied_cells, occupancy.size());
        result.max_occupied_cells = std::max(result.max_occupied_cells, occupancy.size());
        const long double cell_mass = 1.0L / static_cast<long double>(occupancy.size());
        for (std::size_t selected = 0; selected < informed.size(); ++selected) {
            accumulated[selected] += cell_mass / static_cast<long double>(occupancy.find(sample_cells[selected])->second);
        }
    }

    // Dividing the sum across origins by its total is equivalent to averaging
    // the individually normalized lattices and then normalizing once more.
    long double total = 0.0L;
    for (std::vector<long double>::const_iterator it_weight=(accumulated).begin();it_weight!=(accumulated).end();++it_weight) for (bool once_weight=true;once_weight;) for (const std::vector<long double>::value_type weight=*it_weight;once_weight;once_weight=false) total += weight;
    if (!without_sgems::finite(total) || total <= 0.0L) {
        throw std::runtime_error("Cell declustering produced an invalid total weight.");
    }
    long double scaled_mean = 0.0L;
    for (std::size_t selected = 0; selected < informed.size(); ++selected) {
        const long double normalized = accumulated[selected] / total;
        result.weights[informed[selected]] = static_cast<double>(normalized);
        if (largest_absolute_value > 0.0) {
            scaled_mean += normalized *
                (static_cast<long double>(data.nodes[informed[selected]].value)
                    / static_cast<long double>(largest_absolute_value));
        }
    }
    // Scaling prevents an intermediate weighted sum overflowing for finite
    // values near the floating-point limit. The mean is a convex combination.
    scaled_mean = std::max(-1.0L, std::min(1.0L, scaled_mean));
    result.weighted_mean = static_cast<double>(scaled_mean * static_cast<long double>(largest_absolute_value));
    return result;
}

}
