#ifndef WITHOUT_SGEMS_CELL_DECLUSTERING_H
#define WITHOUT_SGEMS_CELL_DECLUSTERING_H

#include "point_set.h"
#include <array>
#include <cstddef>
#include <vector>

namespace without_sgems {

struct CellDeclusteringOptions {
    std::array<double, 3> cell_size{{1.0, 1.0, 1.0}};
    std::array<double, 3> origin{{0.0, 0.0, 0.0}};
    bool origin_at_data_min = true;
    int shifts_per_axis = 3;
    // Zero follows the input dimensions; 2 explicitly groups a 3-D well set in XY.
    int dimensions = 0;
    // Optional reproducible fractional origins in [0,1), replacing regular shifts.
    std::vector<std::array<double, 3>> origin_offsets;
};

struct CellDeclusteringResult {
    // Source-row aligned; uninformed rows have zero weight. Valid weights sum to one.
    std::vector<double> weights;
    std::array<double, 3> origin{{0.0, 0.0, 0.0}};
    std::size_t valid_count = 0;
    std::size_t origin_count = 0;
    std::size_t min_occupied_cells = 0;
    std::size_t max_occupied_cells = 0;
    double weighted_mean = 0.0;
};

// Independent implementation of cell declustering. At each lattice origin an
// occupied cell receives equal total mass, shared equally among its informed
// samples. Origins use every combination of j / shifts_per_axis of a cell
// width along the active axes (j = 0, ..., shifts_per_axis - 1). The origin
// count is limited to 4096. Explicit fractional origin_offsets replace the
// regular schedule; dimensions may restrict grouping to XY for 3-D wells.
// Cells are half-open [lower, upper).
CellDeclusteringResult calculate_cell_declustering(const PointSetData& data,
    const CellDeclusteringOptions& options);

}
#endif
