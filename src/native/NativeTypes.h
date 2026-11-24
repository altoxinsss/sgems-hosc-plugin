#ifndef HOSC_NATIVE_TYPES_H
#define HOSC_NATIVE_TYPES_H
#include <algorithm>
#include <cmath>
#include <float.h>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace without_sgems {
template <class T> struct Triple {
  T v[3];
  Triple(T x = T(), T y = T(), T z = T()) {
    v[0] = x;
    v[1] = y;
    v[2] = z;
  }
  T &operator[](size_t i) { return v[i]; }
  const T &operator[](size_t i) const { return v[i]; }
  T *begin() { return v; }
  T *end() { return v + 3; }
  const T *begin() const { return v; }
  const T *end() const { return v + 3; }
  bool operator<(const Triple &b) const {
    for (int i = 0; i < 3; ++i) {
      if (v[i] < b[i])
        return true;
      if (b[i] < v[i])
        return false;
    }
    return false;
  }
};
typedef Triple<double> Vec3;
inline bool finite(double x) { return _finite(x) != 0; }
inline double nan() { return std::numeric_limits<double>::quiet_NaN(); }
inline double hypot3(double x, double y, double z) { return _hypot(_hypot(x, y), z); }
inline long double clamp_mean(long double x, long double a, long double b) {
  return std::max(a, std::min(b, x));
}
struct PointSample {
  int index;
  double x, y, z, value, external_weight;
  bool informed;
  PointSample() : index(0), x(0), y(0), z(0), value(nan()), external_weight(nan()), informed(false) {}
};
struct PointSetData {
  std::string title, value_label;
  int dimensions;
  std::vector<PointSample> nodes;
  PointSetData() : dimensions(3) {}
};
struct CellDeclusteringOptions {
  Vec3 cell_size, origin;
  bool origin_at_data_min;
  int shifts_per_axis, dimensions;
  std::vector<Vec3> origin_offsets;
  CellDeclusteringOptions()
      : cell_size(1, 1, 1), origin_at_data_min(true), shifts_per_axis(3), dimensions(0) {}
};
struct CellDeclusteringResult {
  std::vector<double> weights;
  Vec3 origin;
  size_t valid_count, origin_count, min_occupied_cells, max_occupied_cells;
  double weighted_mean;
  CellDeclusteringResult()
      : valid_count(0), origin_count(0), min_occupied_cells(0), max_occupied_cells(0),
        weighted_mean(0) {}
};
struct DeclusteringDomain {
  int dimensions;
  Vec3 lower, upper;
  DeclusteringDomain() : dimensions(0) {}
};
struct KrigingDeclusteringOptions {
  DeclusteringDomain domain;
  std::string covariance_model;
  Vec3 ranges;
  double nugget_fraction;
  Triple<int> integration_cells;
  KrigingDeclusteringOptions() : nugget_fraction(0), integration_cells(0, 0, 1) {}
};
struct SpatialDeclusteringResult {
  std::vector<double> weights;
  size_t valid_count, unique_site_count, integration_point_count;
  double weighted_mean, partition_relative_error, kriging_lagrange_multiplier, kriging_max_residual;
  double kriging_min_cholesky_pivot, min_site_weight, max_site_weight;
  SpatialDeclusteringResult()
      : valid_count(0), unique_site_count(0), integration_point_count(0), weighted_mean(0),
        partition_relative_error(0), kriging_lagrange_multiplier(0), kriging_max_residual(0),
        kriging_min_cholesky_pivot(0), min_site_weight(0), max_site_weight(0) {}
};
CellDeclusteringResult calculate_cell_declustering(const PointSetData &,
                                                   const CellDeclusteringOptions &);
SpatialDeclusteringResult calculate_polygonal_declustering(const PointSetData &,
                                                           const DeclusteringDomain &);
SpatialDeclusteringResult calculate_kriging_declustering(const PointSetData &,
                                                         const KrigingDeclusteringOptions &);

class PointSpatialIndex {
public:
  explicit PointSpatialIndex(const PointSetData &);
  std::vector<size_t> query_box(const Vec3 &, const Vec3 &) const;
  std::vector<size_t> query_radius(const Vec3 &, double) const;

private:
  struct Branch {
    size_t point;
    int axis, left, right;
    Branch(size_t p, int a) : point(p), axis(a), left(-1), right(-1) {}
  };
  const PointSetData &data_;
  std::vector<Branch> branches_;
  int root_;
  int build(std::vector<size_t> &, size_t, size_t, int);
  void query(int, const Vec3 &, const Vec3 &, std::vector<size_t> &) const;
};
} // namespace without_sgems
#endif
