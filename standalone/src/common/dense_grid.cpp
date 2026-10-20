#include "dense_grid.h"

#include <stdexcept>

namespace without_sgems {

DenseGrid build_dense_grid(const GridData& source) {
    DenseGrid dense;
    dense.nx = source.nx;
    dense.ny = source.ny;
    dense.nz = source.nz;
    if (dense.nx <= 0 || dense.ny <= 0 || dense.nz <= 0) {
        throw std::runtime_error("Cannot densify a grid with a non-positive dimension.");
    }
    const size_t total = static_cast<size_t>(dense.nx) * dense.ny * dense.nz;
    dense.values.assign(total, 0.0);
    dense.present.assign(total, 0);
    for (size_t i = 0; i < source.nodes.size(); ++i) {
        const GridNode& node = source.nodes[i];
        if (node.x < 0 || node.x >= dense.nx || node.y < 0 || node.y >= dense.ny ||
            node.z < 0 || node.z >= dense.nz) {
            throw std::runtime_error("Encountered an out-of-range node while densifying the source grid.");
        }
        const size_t linear = static_cast<size_t>(node.z) * dense.nx * dense.ny +
            static_cast<size_t>(node.y) * dense.nx + static_cast<size_t>(node.x);
        dense.values[linear] = node.value;
        dense.present[linear] = 1;
    }
    return dense;
}

bool DenseGridView::contains(int x, int y, int z) const {
    return x >= 0 && x < grid_.nx && y >= 0 && y < grid_.ny && z >= 0 && z < grid_.nz;
}

size_t DenseGridView::index(int x, int y, int z) const {
    if (!contains(x, y, z)) {
        throw std::out_of_range("DenseGridView coordinate is outside the source grid.");
    }
    return static_cast<size_t>(z) * grid_.nx * grid_.ny +
        static_cast<size_t>(y) * grid_.nx + static_cast<size_t>(x);
}

bool DenseGridView::present(int x, int y, int z) const {
    return contains(x, y, z) && grid_.present[index(x, y, z)] != 0;
}

double DenseGridView::value(int x, int y, int z) const {
    return grid_.values[index(x, y, z)];
}

} // namespace without_sgems
