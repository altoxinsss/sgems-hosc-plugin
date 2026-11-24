#ifndef WITHOUT_SGEMS_POINT_SET_H
#define WITHOUT_SGEMS_POINT_SET_H

#include "config.h"
#include <array>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

namespace without_sgems {

struct PointSample {
    int index = 0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double value = std::numeric_limits<double>::quiet_NaN();
    bool informed = false;
    double external_weight = std::numeric_limits<double>::quiet_NaN();
};

struct PointSetData {
    std::string title;
    std::string value_label;
    int dimensions = 3;
    std::vector<PointSample> nodes;
};

// Point-set coordinates are physical doubles: no shifting, gridding or rounding.
bool uses_point_set(const IniConfig& config);
PointSetData load_point_set(const std::filesystem::path& path,
    const std::string& source_property = std::string(),
    int coordinate_dimensions = 0,
    double missing_value = std::numeric_limits<double>::quiet_NaN(),
    const std::string& weight_property = std::string());
void write_point_set(const std::filesystem::path& path, const std::string& title,
    const std::string& value_label, const std::vector<PointSample>& nodes);

// Queries return source-vector indices in input order, with informed values only.
// Keep the input data alive and unchanged for the lifetime of this index.
class PointSpatialIndex {
public:
    explicit PointSpatialIndex(const PointSetData& data);
    std::vector<size_t> query_box(const std::array<double, 3>& lower,
        const std::array<double, 3>& upper) const;
    std::vector<size_t> query_radius(const std::array<double, 3>& center,
        double radius) const;
private:
    struct Branch { size_t point; int axis; int left; int right; };
    const PointSetData& data_;
    std::vector<Branch> branches_;
    int root_ = -1;
    int build(std::vector<size_t>& ids, size_t first, size_t last, int depth);
    void query(int branch, const std::array<double, 3>& lower,
        const std::array<double, 3>& upper, std::vector<size_t>& result) const;
};

}
#endif
