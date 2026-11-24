#include "point_set.h"
#include "external_weight_io.h"
#include "text.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>

namespace without_sgems {
namespace {
std::vector<std::string> tokens(std::string line) {
    for (char& c : line) if (c == ',') c = ' ';
    std::istringstream stream(line);
    std::vector<std::string> result;
    std::string token;
    while (stream >> token) result.push_back(token);
    return result;
}

double number(const std::string& text, const std::string& context) {
    try {
        size_t used = 0;
        const double value = std::stod(text, &used);
        if (used != text.size()) throw std::invalid_argument("trailing text");
        return value;
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid numeric value '" + text + "' in " + context);
    }
}

bool integer_token(const std::string& token, int& value) {
    try {
        size_t used = 0;
        const long long parsed = std::stoll(token, &used);
        if (used != token.size() || parsed < 0 || parsed > std::numeric_limits<int>::max()) return false;
        value = static_cast<int>(parsed);
        return true;
    } catch (const std::exception&) { return false; }
}

double coordinate(const PointSample& p, int axis) {
    return axis == 0 ? p.x : axis == 1 ? p.y : p.z;
}
}

bool uses_point_set(const IniConfig& config) {
    bool specified = false;
    bool pointset = false;
    for (const std::string key : {"source_type", "grid_type", "data_type"}) {
        if (!config.has("general", key)) continue;
        const std::string value = to_lower(trim(config.get("general", key)));
        bool next;
        if (value == "point_set" || value == "pointset" || value == "irregular") next = true;
        else if (value == "cartesian" || value == "regular" || value == "grid") next = false;
        else throw std::runtime_error(key + " must be point_set or cartesian.");
        if (specified && pointset != next) throw std::runtime_error("Conflicting input data types in [general].");
        pointset = next;
        specified = true;
    }
    if (config.has("general", "reg_irreg")) {
        const std::string value = trim(config.get("general", "reg_irreg"));
        if (value != "0" && value != "1") throw std::runtime_error("reg_irreg must be 0 (point set) or 1 (Cartesian).");
        const bool next = value == "0";
        if (specified && pointset != next) throw std::runtime_error("source_type and reg_irreg disagree.");
        pointset = next;
    }
    return pointset;
}

PointSetData load_point_set(const std::filesystem::path& path,
    const std::string& source_property, int coordinate_dimensions, double missing_value,
    const std::string& weight_property) {
    if (coordinate_dimensions < 0 || coordinate_dimensions > 3)
        throw std::runtime_error("coordinate_dimensions must be 0, 1, 2 or 3.");
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Unable to open point set: " + path.string());
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(input, line)) lines.push_back(trim(line));
    if (lines.size() < 2) throw std::runtime_error("Point-set file is too short: " + path.string());
    if (lines[0].compare(0, 3, "\xEF\xBB\xBF") == 0) lines[0].erase(0, 3);

    PointSetData data;
    data.title = lines[0];
    const auto first = tokens(lines[0]);
    const bool named_header = first.size() >= 2 &&
        (to_lower(first[0]) == "x" || to_lower(first[0]) == "easting" || to_lower(first[0]) == "longitude");
    const auto second = tokens(lines[1]);
    int first_count = 0, second_count = 0;
    const bool legacy = !named_header && second.size() == 2 && integer_token(second[0], first_count) &&
        integer_token(second[1], second_count);
    const bool gslib = second.size() == 1 && integer_token(second[0], first_count);
    size_t start = 0;
    int value_column = 0, columns = 0;
    int weight_column = -1;
    int expected_count = -1;
    int supplied_coordinates = 0;
    if (legacy) {
        if (!weight_property.empty()) throw std::runtime_error("External weights require named columns (GSLIB or CSV), not legacy HOSC input.");
        if (first_count <= 0 || lines.size() < 5) throw std::runtime_error("Invalid legacy point-set header.");
        expected_count = first_count;
        data.dimensions = coordinate_dimensions == 0 ? 3 : coordinate_dimensions;
        supplied_coordinates = 3;
        data.value_label = source_property.empty() ? "value" : source_property;
        start = 4;
        columns = 5;
        value_column = 4;
    } else {
        std::vector<std::string> labels;
        if (gslib) {
            columns = first_count;
            if (columns < 2 || lines.size() <= 2 + static_cast<size_t>(columns))
                throw std::runtime_error("Invalid GSLIB point-set header or missing data.");
            for (int i = 0; i < columns; ++i) labels.push_back(lines[2 + static_cast<size_t>(i)]);
            start = 2 + static_cast<size_t>(columns);
        } else {
            labels = tokens(lines[0]);
            columns = static_cast<int>(labels.size());
            start = 1;
            if (columns < 2) throw std::runtime_error("Expected GSLIB, legacy HOSC or a named-column point set.");
        }
        int inferred = 0;
        const char* axes[] = {"x", "y", "z"};
        while (inferred < 3 && inferred < columns - 1 && to_lower(labels[inferred]) == axes[inferred]) ++inferred;
        const int named_coordinates = inferred;
        if (inferred == 0) inferred = std::min(3, columns - 1);
        data.dimensions = coordinate_dimensions == 0 ? inferred : coordinate_dimensions;
        supplied_coordinates = std::max(data.dimensions, named_coordinates);
        if (data.dimensions >= columns) throw std::runtime_error("coordinate_dimensions leaves no property column.");
        value_column = columns - 1;
        if (!source_property.empty()) {
            value_column = -1;
            for (int i = supplied_coordinates; i < columns; ++i)
                if (to_lower(labels[i]) == to_lower(source_property)) {
                    if (value_column >= 0) throw std::runtime_error("Duplicate source_property column name.");
                    value_column = i;
                }
            if (value_column < 0) throw std::runtime_error("source_property '" + source_property + "' not found in " + path.string());
        }
        data.value_label = labels[value_column];
        weight_column = external_weight_column(labels, weight_property, supplied_coordinates);
    }
    for (size_t row = start; row < lines.size(); ++row) {
        if (lines[row].empty() || lines[row][0] == '#') continue;
        const auto parts = tokens(lines[row]);
        if (parts.size() != static_cast<size_t>(columns))
            throw std::runtime_error("Wrong column count at point-set line " + std::to_string(row + 1));
        if (data.nodes.size() >= static_cast<size_t>(std::numeric_limits<int>::max()))
            throw std::runtime_error("Too many point-set rows.");
        PointSample p;
        p.index = static_cast<int>(data.nodes.size()) + 1;
        const std::string context = path.string() + ":" + std::to_string(row + 1);
        const int offset = legacy ? 1 : 0;
        if (legacy) {
            const double id = number(parts[0], context);
            if (!std::isfinite(id) || id < std::numeric_limits<int>::min() || id > std::numeric_limits<int>::max() || std::floor(id) != id)
                throw std::runtime_error("Invalid point ID in " + context);
            p.index = static_cast<int>(id);
        }
        p.x = number(parts[offset], context);
        p.y = supplied_coordinates >= 2 ? number(parts[offset + 1], context) : 0.0;
        p.z = supplied_coordinates >= 3 ? number(parts[offset + 2], context) : 0.0;
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
            throw std::runtime_error("Point coordinates must be finite in " + context);
        if (data.dimensions < 3 && p.z != 0.0)
            throw std::runtime_error("Nonzero z coordinate conflicts with coordinate_dimensions in " + context);
        if (data.dimensions < 2 && p.y != 0.0)
            throw std::runtime_error("Nonzero y coordinate conflicts with coordinate_dimensions in " + context);
        p.value = number(parts[value_column], context);
        p.informed = std::isfinite(p.value) && (!std::isfinite(missing_value) || p.value != missing_value);
        if (!p.informed) p.value = std::numeric_limits<double>::quiet_NaN();
        if (weight_column >= 0) p.external_weight = external_weight_token(parts[weight_column]);
        data.nodes.push_back(p);
    }
    if (data.nodes.empty()) throw std::runtime_error("Point set contains no data rows.");
    if (expected_count >= 0 && data.nodes.size() != static_cast<size_t>(expected_count))
        throw std::runtime_error("Legacy point-set row count disagrees with ND.");
    return data;
}

void write_point_set(const std::filesystem::path& path, const std::string& title,
    const std::string& value_label, const std::vector<PointSample>& nodes) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("Unable to create point-set output: " + path.string());
    output << title << "\n4\nx\ny\nz\n" << value_label << "\n" << std::setprecision(17);
    for (const auto& p : nodes) {
        output << p.x << '\t' << p.y << '\t' << p.z << '\t';
        if (p.informed && std::isfinite(p.value)) output << p.value;
        else output << "nan";
        output << '\n';
    }
    if (!output) throw std::runtime_error("Failed to write point-set output: " + path.string());
}

PointSpatialIndex::PointSpatialIndex(const PointSetData& data) : data_(data) {
    if (data.nodes.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::runtime_error("Too many points for spatial index.");
    std::vector<size_t> ids;
    for (size_t i = 0; i < data.nodes.size(); ++i) {
        const auto& p = data.nodes[i];
        if (!std::isfinite(p.x) || !std::isfinite(p.y) || !std::isfinite(p.z))
            throw std::runtime_error("Point coordinates must be finite.");
        if (p.informed && std::isfinite(p.value)) ids.push_back(i);
    }
    branches_.reserve(ids.size());
    root_ = build(ids, 0, ids.size(), 0);
}

int PointSpatialIndex::build(std::vector<size_t>& ids, size_t first, size_t last, int depth) {
    if (first == last) return -1;
    const size_t mid = first + (last - first) / 2;
    const int axis = depth % 3;
    std::nth_element(ids.begin() + first, ids.begin() + mid, ids.begin() + last,
        [&](size_t a, size_t b) {
            const double va = coordinate(data_.nodes[a], axis), vb = coordinate(data_.nodes[b], axis);
            return va < vb || (va == vb && a < b);
        });
    const int current = static_cast<int>(branches_.size());
    branches_.push_back({ids[mid], axis, -1, -1});
    const int left = build(ids, first, mid, depth + 1);
    const int right = build(ids, mid + 1, last, depth + 1);
    branches_[current].left = left;
    branches_[current].right = right;
    return current;
}

void PointSpatialIndex::query(int branch, const std::array<double, 3>& lower,
    const std::array<double, 3>& upper, std::vector<size_t>& result) const {
    if (branch < 0) return;
    const auto& b = branches_[branch];
    const auto& p = data_.nodes[b.point];
    if (p.x >= lower[0] && p.x <= upper[0] && p.y >= lower[1] && p.y <= upper[1] && p.z >= lower[2] && p.z <= upper[2])
        result.push_back(b.point);
    const double value = coordinate(p, b.axis);
    if (lower[b.axis] <= value) query(b.left, lower, upper, result);
    if (upper[b.axis] >= value) query(b.right, lower, upper, result);
}

std::vector<size_t> PointSpatialIndex::query_box(const std::array<double, 3>& lower,
    const std::array<double, 3>& upper) const {
    for (size_t i = 0; i < 3; ++i)
        if (std::isnan(lower[i]) || std::isnan(upper[i]) || lower[i] > upper[i])
            throw std::runtime_error("Invalid spatial query bounds.");
    std::vector<size_t> result;
    query(root_, lower, upper, result);
    std::sort(result.begin(), result.end());
    return result;
}

std::vector<size_t> PointSpatialIndex::query_radius(const std::array<double, 3>& center, double radius) const {
    if (!std::isfinite(radius) || radius < 0.0) throw std::runtime_error("Spatial query radius must be finite and nonnegative.");
    std::array<double, 3> lower, upper;
    for (size_t i = 0; i < 3; ++i) {
        if (!std::isfinite(center[i])) throw std::runtime_error("Spatial query center must be finite.");
        lower[i] = center[i] - radius;
        upper[i] = center[i] + radius;
    }
    auto result = query_box(lower, upper);
    result.erase(std::remove_if(result.begin(), result.end(), [&](size_t i) {
        const auto& p = data_.nodes[i];
        return std::hypot(p.x - center[0], p.y - center[1], p.z - center[2]) > radius;
    }), result.end());
    return result;
}
}
