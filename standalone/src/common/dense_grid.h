#ifndef WITHOUT_SGEMS_DENSE_GRID_H
#define WITHOUT_SGEMS_DENSE_GRID_H

#include "grid_io.h"

#include <vector>

namespace without_sgems {

struct DenseGrid {
    int nx = 0;
    int ny = 0;
    int nz = 0;
    std::vector<double> values;
    std::vector<unsigned char> present;
};

struct GridBounds {
    int x0 = 0;
    int y0 = 0;
    int z0 = 0;
    int nx = 0;
    int ny = 0;
    int nz = 0;
};

DenseGrid build_dense_grid(const GridData& source);

class DenseGridView {
public:
    explicit DenseGridView(const DenseGrid& grid) : grid_(grid) {}

    bool contains(int x, int y, int z) const;
    bool present(int x, int y, int z) const;
    double value(int x, int y, int z) const;

    int nx() const { return grid_.nx; }
    int ny() const { return grid_.ny; }
    int nz() const { return grid_.nz; }

private:
    const DenseGrid& grid_;
    size_t index(int x, int y, int z) const;
};

} // namespace without_sgems

#endif
