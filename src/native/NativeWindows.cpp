#include "NativeInternal.h"
#include "WindowGeometry.h"

namespace hosc_native {
namespace {
struct Ellipsoid : hosc_window_geometry::Shape {
  Ellipsoid(const Config& c,int dimensions)
    : hosc_window_geometry::Shape(c.get("sliding_window","neighborhood_shape","ellipsoid"),
        c.list("sliding_window","neighborhood_radii"),
        c.list("sliding_window","neighborhood_angles","90,0,0"),dimensions) {}
};
struct Result {
  Vec3 coord;
  double value;
  size_t count;
  long double weight_sum, weight_squares;
  bool complete;
  Result() : value(without_sgems::nan()), count(0), weight_sum(0), weight_squares(0), complete(true) {}
  void add_weight(double w) { weight_sum += w; weight_squares += (long double)w * w; }
};
double percentile(const std::vector<double> &sorted, double pct) {
  double pos = pct * .01 * (sorted.size() - 1);
  size_t lo = (size_t)floor(pos), hi = (size_t)ceil(pos);
  return sorted[lo] + (pos - lo) * (sorted[hi] - sorted[lo]);
}
} // namespace
void calculate_windows(const Config &c, const Source &s, const std::string &dir) {
  const int order = integer(number(c.get("sliding_window", "moment_cumulant_order")), 1);
  bool moment = lower(c.get("general", "calculate_target")) == "moment",
       bell = lower(c.get("general", "calculate_algorithm")).find("bell") != std::string::npos;
  bool invert = c.flag("sliding_window", "invert_window_output"),
       clip = c.flag("sliding_window", "enable_percentile_clipping");
  std::vector<Result> results;
  hosc_window_geometry::CenterGrid centers;
  std::vector<int> powers(1, order);
  const bool weighted = c.get("weighted_moment", "weight_mode", "uniform") != "uniform";
  std::ostringstream weight_meta;
  weight_meta << std::setprecision(17);
  // Shared with templates: compute source-point weights once, then select
  // and normalize them on each window's informed sample support.
  std::vector<double> source_weights = point_weights(c, s.data, weight_meta);
  if (s.points) {
    Ellipsoid e(c, s.data.dimensions);
    without_sgems::PointSpatialIndex index(s.data);
    bool includes = c.flag("sliding_window", "includes_center", true);
    size_t minimum = integer(number(c.get("sliding_window", "min_samples", "1")), 1);
    const std::string mode=c.get("sliding_window","center_mode","sample_points");
    const std::string bounds=c.get("sliding_window","grid_bounds_mode","auto");
    const std::vector<double> empty;
    centers.configure(s.data,e,mode,bounds,
        mode=="regular_grid"&&bounds=="manual"?c.list("sliding_window","grid_bounds_min"):empty,
        mode=="regular_grid"&&bounds=="manual"?c.list("sliding_window","grid_bounds_max"):empty,
        mode=="regular_grid"?c.list("sliding_window","grid_start"):empty,
        mode=="regular_grid"?c.list("sliding_window","grid_spacing"):empty);
    if(centers.regular&&!includes)
        throw std::runtime_error("includes_center=false applies only to sample_points mode.");
    for (size_t i = 0; i < centers.count; ++i) {
      const PointSample center = centers.regular?centers.point<PointSample>(i):s.data.nodes[i];
      Result r;
      r.coord = Vec3(center.x, center.y, center.z);
      Vec3 lo, hi;
      for (int k = 0; k < 3; ++k) {
        double pad = 1e-12 * std::max(1.0, e.extent[k]);
        lo[k] = r.coord[k] - e.extent[k] - pad;
        hi[k] = r.coord[k] + e.extent[k] + pad;
      }
      std::vector<size_t> ids = index.query_box(lo, hi);
      Values samples;
      std::vector<double> weights;
      for (size_t j = 0; j < ids.size(); ++j)
        if ((centers.regular || includes || ids[j] != i) && e.contains(s.data.nodes[ids[j]], center)) {
          samples.push_back(std::vector<double>(1, s.data.nodes[ids[j]].value));
          r.add_weight(source_weights[ids[j]]);
          if (weighted) weights.push_back(source_weights[ids[j]]);
        }
      r.count = samples.size();
      if (r.count >= minimum)
        r.value = evaluate(samples, weights, powers, moment, bell);
      results.push_back(r);
    }
  } else {
    Vec3 window = vector3(c.list("sliding_window", "window_size")),
         step = vector3(c.list("sliding_window", "move_step"));
    int dims[] = {s.nx, s.ny, s.nz}, w[3], st[3], outdim[3];
    for (int k = 0; k < 3; ++k) {
      w[k] = integer(window[k], 1, dims[k]);
      st[k] = integer(step[k], k == 2 && s.nz == 1 ? 0 : 1);
      if (k == 2 && s.nz == 1 && (w[k] != 1 || st[k] != 0))
        throw std::runtime_error("2D window Z must be 1 and step Z must be 0.");
      outdim[k] = st[k] ? 1 + (dims[k] - w[k]) / st[k] : 1;
    }
    for (int az = 0; az < outdim[2]; ++az)
      for (int ay = 0; ay < outdim[1]; ++ay)
        for (int ax = 0; ax < outdim[0]; ++ax) {
          int x0 = ax * st[0], y0 = ay * st[1], z0 = az * st[2];
          Result r;
          r.coord = Vec3(x0 + (w[0] - 1) / 2, y0 + (w[1] - 1) / 2, z0 + (w[2] - 1) / 2);
          Values samples;
          std::vector<double> weights;
          for (int z = 0; z < w[2]; ++z)
            for (int y = 0; y < w[1]; ++y)
              for (int x = 0; x < w[0]; ++x) {
                const size_t id = ((size_t)(z0 + z) * s.ny + y0 + y) * s.nx + x0 + x;
                const PointSample &p = s.data.nodes[id];
                if (p.informed) {
                  samples.push_back(std::vector<double>(1, p.value));
                  r.add_weight(source_weights[id]);
                  if (weighted) weights.push_back(source_weights[id]);
                }
                else
                  r.complete = false;
              }
          r.count = samples.size();
          if (r.complete)
            r.value = evaluate(samples, weights, powers, moment, bell);
          results.push_back(r);
        }
  }
  std::vector<double> finite_values;
  for (size_t i = 0; i < results.size(); ++i)
    if (without_sgems::finite(results[i].value)) {
      if (invert)
        results[i].value = -results[i].value;
      finite_values.push_back(results[i].value);
    }
  if (clip && !finite_values.empty()) {
    double low = number(c.get("sliding_window", "percentile_vmin", "2")),
           high = number(c.get("sliding_window", "percentile_vmax", "98"));
    if (low > high)
      std::swap(low, high);
    if (s.points && (low < 0 || high > 100))
      throw std::runtime_error("Clipping percentiles must be in [0,100].");
    low = std::max(0.0, std::min(100.0, low));
    high = std::max(0.0, std::min(100.0, high));
    std::sort(finite_values.begin(), finite_values.end());
    double lo = percentile(finite_values, low), hi = percentile(finite_values, high);
    for (size_t i = 0; i < results.size(); ++i)
      if (without_sgems::finite(results[i].value))
        results[i].value = std::max(lo, std::min(hi, results[i].value));
  }
  std::ofstream out((dir + "\\result.txt").c_str()),
      support((dir + "\\result.window_support.csv").c_str()),
      meta((dir + "\\result.meta.txt").c_str());
  out << "HOSC window result\n4\nx\ny\nz\nresult\n" << std::setprecision(17);
  support << "x,y,z,weight_sum,weight_squared_sum,effective_sample_size,sample_count,map_included,"
             "valid,valid_reason\n"
          << std::setprecision(17);
  for (size_t i = 0; i < results.size(); ++i) {
    const Result &r = results[i];
    bool valid = without_sgems::finite(r.value);
    out << r.coord[0] << ' ' << r.coord[1] << ' ' << r.coord[2] << ' ';
    write_value(out, r.value);
    out << '\n';
    long double sum = s.points || valid ? r.weight_sum : 0;
    long double squares = s.points || valid ? r.weight_squares : 0;
    double eff = squares > 0 ? (double)(sum * sum / squares) : 0;
    support << r.coord[0] << ',' << r.coord[1] << ',' << r.coord[2] << ',' << (double)sum << ',' << (double)squares
            << ',' << eff << ',' << r.count << ',' << ((s.points || valid) ? 1 : 0) << ','
            << (valid ? 1 : 0) << ','
            << (valid          ? "valid"
                : !r.complete  ? "incomplete_window"
                : r.count == 0 ? "empty_neighborhood"
                               : "minimum_samples")
            << '\n';
  }
  if (weighted) {
    std::ofstream pointout((dir + "\\result.point_weights.csv").c_str());
    pointout << "node_id,x,y,z,value,informed,weight\n" << std::setprecision(17);
    for (size_t i = 0; i < s.data.nodes.size(); ++i) {
      const PointSample &p = s.data.nodes[i];
      pointout << p.index << ',' << p.x << ',' << p.y << ',' << p.z << ',';
      write_value(pointout, p.value);
      pointout << ',' << (p.informed ? 1 : 0) << ',' << source_weights[i] << '\n';
    }
    if (!pointout) throw std::runtime_error("Cannot write window source point weights.");
  }
  meta << weight_meta.str()
       << "normalized_weights=true\npoint_weight_scope=source_dataset\nwindow_weight_normalization=selected_informed_samples\n"
       << "native_vc71_in_process=true\norder=" << order << "\nevaluator="
       << (moment ? "raw_moment"
           : bell ? "univariate_recurrence"
                  : "AVL")
       << "\ncumulant_centering=within_complete_support\n";
  if(s.points) {
    centers.metadata(meta);
    meta<<"neighborhood_shape="<<c.get("sliding_window","neighborhood_shape","ellipsoid")<<'\n'
        <<"neighborhood_radii="<<c.get("sliding_window","neighborhood_radii")<<'\n'
        <<"neighborhood_angles="<<c.get("sliding_window","neighborhood_angles","90,0,0")<<'\n';
  }
  if (!out || !support || !meta)
    throw std::runtime_error("Cannot write moving-window results.");
}
} // namespace hosc_native
