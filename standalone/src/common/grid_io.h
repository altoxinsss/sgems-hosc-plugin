#ifndef WITHOUT_SGEMS_GRID_IO_H
#define WITHOUT_SGEMS_GRID_IO_H

#include <filesystem>
#include <string>
#include <vector>

namespace without_sgems {

enum class GridTextFormat {
    LegacyHosc,
    ColumnarNamed,
};

struct GridNode {
    int index;
    int x;
    int y;
    int z;
    double value;
};

struct GridData {
    int nd;
    int nvar;
    int nx;
    int ny;
    int nz;
    int raw_min_x;
    int raw_min_y;
    int raw_min_z;
    GridTextFormat text_format;
    std::string title;
    std::string x_label;
    std::string y_label;
    std::string z_label;
    std::string value_label;
    std::vector<GridNode> nodes;
};

// For a columnar grid, source_property selects a named value column.  Set
// coordinate_dimensions for layouts whose coordinate columns cannot be
// inferred from their column count (for example x,y,porosity,permeability).
GridData load_hosc_grid(
    const std::filesystem::path& path,
    const std::string& source_property = std::string(),
    int coordinate_dimensions = 0);
double** create_vendor_grid_array(const GridData& grid);
void free_vendor_matrix(double** matrix, int rows);
void write_text_grid(
    const std::filesystem::path& path,
    const GridData& style,
    const std::string& title,
    const std::string& value_label,
    const std::vector<GridNode>& nodes);

}

#endif
