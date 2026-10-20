#include "grid_io.h"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <stdexcept>
#include <vector>

namespace without_sgems {

namespace {

std::vector<std::string> read_all_lines(const std::filesystem::path& path) {
    std::ifstream input(path.c_str(), std::ios::in);
    if (!input) {
        throw std::runtime_error("Unable to open source grid file: " + path.string());
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line)) {
        lines.push_back(line);
    }
    return lines;
}

std::vector<std::string> split_whitespace(const std::string& line) {
    std::vector<std::string> tokens;
    std::istringstream stream(line);
    std::string token;
    while (stream >> token) {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<int> parse_int_tokens(const std::string& line) {
    std::vector<std::string> tokens = split_whitespace(line);
    std::vector<int> values;
    values.reserve(tokens.size());

    for (size_t i = 0; i < tokens.size(); ++i) {
        char* end = NULL;
        long parsed = std::strtol(tokens[i].c_str(), &end, 10);
        if (end == NULL || *end != '\0') {
            values.clear();
            return values;
        }
        values.push_back(static_cast<int>(parsed));
    }
    return values;
}

std::string lower_ascii(const std::string& value) {
    std::string normalized = value;
    for (size_t i = 0; i < normalized.size(); ++i) {
        normalized[i] = static_cast<char>(std::tolower(static_cast<unsigned char>(normalized[i])));
    }
    return normalized;
}

GridData load_legacy_hosc_grid(
    const std::filesystem::path& path,
    const std::vector<std::string>& lines,
    const std::string& source_property) {
    if (lines.size() < 5) {
        throw std::runtime_error("Invalid legacy HOSC grid header in: " + path.string());
    }

    const std::vector<int> header = parse_int_tokens(lines[1]);
    if (header.size() < 2) {
        throw std::runtime_error("Invalid legacy HOSC grid counts in: " + path.string());
    }

    GridData grid;
    grid.nd = header[0];
    grid.nvar = header[1];
    grid.text_format = GridTextFormat::LegacyHosc;
    grid.title = lines[0];
    grid.x_label = "x";
    grid.y_label = "y";
    grid.z_label = "z";
    grid.value_label = source_property.empty() ? "value" : source_property;
    grid.nodes.reserve(static_cast<size_t>(grid.nd));

    const size_t data_start = 4;
    if (lines.size() < data_start + static_cast<size_t>(grid.nd)) {
        throw std::runtime_error("Legacy HOSC grid row count is shorter than ND in: " + path.string());
    }

    for (int i = 0; i < grid.nd; ++i) {
        const std::vector<std::string> tokens = split_whitespace(lines[data_start + static_cast<size_t>(i)]);
        if (tokens.size() < 5) {
            throw std::runtime_error("Failed to parse legacy HOSC row " + std::to_string(i) + " from: " + path.string());
        }

        GridNode node;
        node.index = std::atoi(tokens[0].c_str());
        node.x = std::atoi(tokens[1].c_str());
        node.y = std::atoi(tokens[2].c_str());
        node.z = std::atoi(tokens[3].c_str());
        node.value = std::atof(tokens[4].c_str());
        grid.nodes.push_back(node);
    }

    return grid;
}

GridData load_columnar_grid(
    const std::filesystem::path& path,
    const std::vector<std::string>& lines,
    const std::string& source_property,
    int coordinate_dimensions) {
    if (lines.size() < 4) {
        throw std::runtime_error("Invalid columnar grid header in: " + path.string());
    }

    const std::vector<int> header = parse_int_tokens(lines[1]);
    if (header.size() != 1) {
        throw std::runtime_error("Invalid column count in: " + path.string());
    }

    const int column_count = header[0];
    if (column_count < 2) {
        throw std::runtime_error("Columnar grid must provide at least one coordinate and one value column: " + path.string());
    }

    const size_t data_start = 2 + static_cast<size_t>(column_count);
    if (lines.size() <= data_start) {
        throw std::runtime_error("Columnar grid contains no data rows: " + path.string());
    }

    GridData grid;
    grid.nvar = 4;
    grid.text_format = GridTextFormat::ColumnarNamed;
    grid.title = lines[0];
    const int inferred_coordinate_dimensions = std::min(3, column_count - 1);
    const int coordinate_count = coordinate_dimensions == 0
        ? inferred_coordinate_dimensions
        : coordinate_dimensions;
    if (coordinate_count < 1 || coordinate_count > 3 || coordinate_count >= column_count) {
        throw std::runtime_error("coordinate_dimensions must leave at least one value column and be between 1 and 3.");
    }
    const size_t header_start = 2;
    grid.x_label = lines[header_start];
    grid.y_label = coordinate_count >= 2 ? lines[header_start + 1] : "y";
    grid.z_label = coordinate_count >= 3 ? lines[header_start + 2] : "z";

    int value_column = column_count - 1;
    if (!source_property.empty()) {
        const std::string requested = lower_ascii(source_property);
        bool found = false;
        for (int column = coordinate_count; column < column_count; ++column) {
            if (lower_ascii(lines[header_start + static_cast<size_t>(column)]) == requested) {
                value_column = column;
                found = true;
                break;
            }
        }
        if (!found) {
            throw std::runtime_error("source_property '" + source_property + "' was not found in: " + path.string());
        }
    }
    grid.value_label = lines[header_start + static_cast<size_t>(value_column)];

    int row_index = 1;
    for (size_t i = data_start; i < lines.size(); ++i) {
        const std::vector<std::string> tokens = split_whitespace(lines[i]);
        if (tokens.empty()) {
            continue;
        }
        if (tokens.size() < static_cast<size_t>(column_count)) {
            throw std::runtime_error("Failed to parse columnar grid row " + std::to_string(i - data_start) + " from: " + path.string());
        }

        GridNode node;
        node.index = row_index++;
        node.x = std::atoi(tokens[0].c_str());
        node.y = coordinate_count >= 2 ? std::atoi(tokens[1].c_str()) : 0;
        node.z = coordinate_count >= 3 ? std::atoi(tokens[2].c_str()) : 0;
        node.value = std::atof(tokens[static_cast<size_t>(value_column)].c_str());
        grid.nodes.push_back(node);
    }

    grid.nd = static_cast<int>(grid.nodes.size());
    if (grid.nd <= 0) {
        throw std::runtime_error("Columnar grid contains no parsed rows: " + path.string());
    }
    return grid;
}

void normalize_coordinates(GridData& grid) {
    if (grid.nodes.empty()) {
        throw std::runtime_error("Grid contains no nodes.");
    }

    int min_x = std::numeric_limits<int>::max();
    int min_y = std::numeric_limits<int>::max();
    int min_z = std::numeric_limits<int>::max();
    int max_x = std::numeric_limits<int>::min();
    int max_y = std::numeric_limits<int>::min();
    int max_z = std::numeric_limits<int>::min();

    for (size_t i = 0; i < grid.nodes.size(); ++i) {
        min_x = std::min(min_x, grid.nodes[i].x);
        min_y = std::min(min_y, grid.nodes[i].y);
        min_z = std::min(min_z, grid.nodes[i].z);
        max_x = std::max(max_x, grid.nodes[i].x);
        max_y = std::max(max_y, grid.nodes[i].y);
        max_z = std::max(max_z, grid.nodes[i].z);
    }

    for (size_t i = 0; i < grid.nodes.size(); ++i) {
        grid.nodes[i].x -= min_x;
        grid.nodes[i].y -= min_y;
        grid.nodes[i].z -= min_z;
    }

    grid.raw_min_x = min_x;
    grid.raw_min_y = min_y;
    grid.raw_min_z = min_z;
    grid.nx = max_x - min_x + 1;
    grid.ny = max_y - min_y + 1;
    grid.nz = max_z - min_z + 1;
}

}

GridData load_hosc_grid(
    const std::filesystem::path& path,
    const std::string& source_property,
    int coordinate_dimensions) {
    const std::vector<std::string> lines = read_all_lines(path);
    if (lines.size() < 2) {
        throw std::runtime_error("Grid file is too short: " + path.string());
    }

    const std::vector<int> second_line = parse_int_tokens(lines[1]);
    GridData grid;
    if (second_line.size() >= 2) {
        grid = load_legacy_hosc_grid(path, lines, source_property);
    }
    else if (second_line.size() == 1) {
        grid = load_columnar_grid(path, lines, source_property, coordinate_dimensions);
    }
    else {
        throw std::runtime_error("Unrecognized grid header format in: " + path.string());
    }

    normalize_coordinates(grid);
    return grid;
}

double** create_vendor_grid_array(const GridData& grid) {
    double** result = new double*[static_cast<size_t>(grid.nd) + 1];
    result[0] = new double[2];
    result[0][0] = static_cast<double>(grid.nd);
    result[0][1] = 4.0;

    for (int i = 0; i < grid.nd; ++i) {
        result[i + 1] = new double[4];
        result[i + 1][0] = static_cast<double>(grid.nodes[static_cast<size_t>(i)].x);
        result[i + 1][1] = static_cast<double>(grid.nodes[static_cast<size_t>(i)].y);
        result[i + 1][2] = static_cast<double>(grid.nodes[static_cast<size_t>(i)].z);
        result[i + 1][3] = grid.nodes[static_cast<size_t>(i)].value;
    }
    return result;
}

void free_vendor_matrix(double** matrix, int rows) {
    if (matrix == NULL) {
        return;
    }

    for (int i = 0; i < rows; ++i) {
        if (matrix[i] != NULL) {
            std::free(matrix[i]);
        }
    }
    std::free(matrix);
}

void write_text_grid(
    const std::filesystem::path& path,
    const GridData& style,
    const std::string& title,
    const std::string& value_label,
    const std::vector<GridNode>& nodes) {
    std::ofstream output(path.c_str(), std::ios::out | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("Unable to create output file: " + path.string());
    }
    output << std::setprecision(17);

    if (style.text_format == GridTextFormat::LegacyHosc) {
        output << title << "\n";
        output << std::string(10, ' ') << nodes.size() << std::string(11, ' ') << 4 << "\n";
        output << " # FOR EACH OF THE NODES :  \"INDEX\" \"X\" \"Y\" \"Z\" \"VAL\"\n";
        output << " # IF THE DATA ARE IN 2D : Z = 0\n";
        for (size_t i = 0; i < nodes.size(); ++i) {
            output << std::string(9, ' ') << nodes[i].index
                   << std::string(11, ' ') << nodes[i].x
                   << std::string(11, ' ') << nodes[i].y
                   << std::string(11, ' ') << nodes[i].z
                   << std::string(11, ' ') << nodes[i].value << "\n";
        }
        return;
    }

    output << title << "\n";
    output << 4 << "\n";
    output << (style.x_label.empty() ? "x" : style.x_label) << "\n";
    output << (style.y_label.empty() ? "y" : style.y_label) << "\n";
    output << (style.z_label.empty() ? "z" : style.z_label) << "\n";
    output << value_label << "\n";
    for (size_t i = 0; i < nodes.size(); ++i) {
        output << nodes[i].x << "\t" << nodes[i].y << "\t" << nodes[i].z << "\t" << nodes[i].value << "\n";
    }
}

}
