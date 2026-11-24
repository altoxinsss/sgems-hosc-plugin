#include "NativeCore.h"
#include "NativeInternal.h"
#include <cctype>
#include <cerrno>
#include <direct.h>
#include <new>

namespace hosc_native {
std::string trim(const std::string &s) {
  size_t a = s.find_first_not_of(" \t\r\n"), b = s.find_last_not_of(" \t\r\n");
  return a == std::string::npos ? "" : s.substr(a, b - a + 1);
}
std::string lower(std::string s) {
  for (size_t i = 0; i < s.size(); ++i)
    s[i] = (char)std::tolower((unsigned char)s[i]);
  return s;
}
double number(const std::string &input) {
  std::string s = trim(input);
  char *end = NULL;
  double v = strtod(s.c_str(), &end);
  if (s.empty() || *end || !without_sgems::finite(v))
    throw std::runtime_error("Expected a finite number: " + input);
  return v;
}
std::vector<double> numbers(const std::string &s) {
  std::vector<double> out;
  if (trim(s).empty())
    return out;
  size_t at = 0;
  while (at <= s.size()) {
    size_t end = s.find(',', at);
    out.push_back(number(s.substr(at, end == std::string::npos ? end : end - at)));
    if (end == std::string::npos)
      break;
    at = end + 1;
  }
  return out;
}
int integer(double v, int minimum, int maximum) {
  if (v != floor(v) || v < minimum || v > maximum)
    throw std::runtime_error("Integer parameter is outside its supported range.");
  return (int)v;
}
Vec3 vector3(const std::vector<double> &v) {
  if (v.size() != 3)
    throw std::runtime_error("Expected three coordinate components.");
  return Vec3(v[0], v[1], v[2]);
}
Config::Config(const std::string &dir, const std::string &filename) : directory(dir) {
  std::ifstream in(path(filename).c_str());
  if (!in)
    throw std::runtime_error("Cannot open calculation configuration.");
  std::string line, section;
  while (std::getline(in, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#' || line[0] == ';')
      continue;
    if (line[0] == '[') {
      size_t end = line.find(']');
      if (end == std::string::npos)
        throw std::runtime_error("Malformed configuration section.");
      section = lower(trim(line.substr(1, end - 1)));
      continue;
    }
    size_t p = line.find('=');
    if (p == std::string::npos)
      throw std::runtime_error("Malformed configuration setting.");
    items[section + "." + lower(trim(line.substr(0, p)))] = trim(line.substr(p + 1));
  }
}
std::string Config::get(const std::string &s, const std::string &k,
                        const std::string &fallback) const {
  std::map<std::string, std::string>::const_iterator i = items.find(lower(s + "." + k));
  return i == items.end() ? fallback : i->second;
}
std::vector<double> Config::list(const std::string &s, const std::string &k,
                                 const std::string &fallback) const {
  return numbers(get(s, k, fallback));
}
bool Config::flag(const std::string &s, const std::string &k, bool fallback) const {
  std::string v = lower(get(s, k, fallback ? "true" : "false"));
  if (v == "true" || v == "1")
    return true;
  if (v == "false" || v == "0")
    return false;
  throw std::runtime_error("Expected a boolean parameter: " + k);
}
std::string Config::path(const std::string &p) const {
  if (p.size() > 1 && (p[1] == ':' || p[0] == '/' || p[0] == '\\'))
    return p;
  return directory + "\\" + p;
}
Source::Source(const Config &c)
    : points(lower(c.get("general", "source_type")) == "point_set"), nx(0), ny(0), nz(0),
      spacing(1, 1, 1) {
  data.dimensions = integer(number(c.get("general", "coordinate_dimensions", "3")), 2, 3);
  origin = vector3(c.list("general", "source_origin", "0,0,0"));
  spacing = vector3(c.list("general", "source_spacing", "1,1,1"));
  std::ifstream in(c.path(c.get("general", "source_grid")).c_str());
  if (!in)
    throw std::runtime_error("Cannot open exported source.");
  std::string line;
  std::getline(in, data.title);
  std::getline(in, line);
  const bool external = c.get("weighted_moment", "weight_mode") == "external";
  if (integer(number(line)) != (external ? 5 : 4))
    throw std::runtime_error("Unexpected exported source column count.");
  for (int i = 0; i < (external ? 5 : 4); ++i)
    std::getline(in, line);
  while (std::getline(in, line)) {
    if (trim(line).empty())
      continue;
    std::istringstream row(line);
    PointSample p;
    std::string val;
    if (!(row >> p.x >> p.y >> p.z >> val))
      throw std::runtime_error("Invalid exported source row.");
    p.value = lower(val) == "nan" ? without_sgems::nan() : number(val);
    p.informed = without_sgems::finite(p.value);
    p.index = (int)data.nodes.size() + 1;
    if (external) {
      if (!(row >> val)) throw std::runtime_error("Missing exported weight column.");
      p.external_weight = lower(val) == "nan" ? without_sgems::nan() : number(val);
    }
    if (!points) {
      nx = std::max(nx, integer(p.x) + 1);
      ny = std::max(ny, integer(p.y) + 1);
      nz = std::max(nz, integer(p.z) + 1);
    }
    data.nodes.push_back(p);
  }
  if (data.nodes.empty())
    throw std::runtime_error("Source has no records.");
  if (points && data.dimensions == 2)
    for (size_t i = 0; i < data.nodes.size(); ++i)
      if (data.nodes[i].z != 0)
        throw std::runtime_error("Nonzero z coordinate conflicts with coordinate_dimensions=2.");
  if (!points) {
    if ((size_t)nx * ny * nz != data.nodes.size())
      throw std::runtime_error("Cartesian source is not a complete lattice.");
    data.dimensions = nz == 1 ? 2 : 3;
    for (size_t i = 0; i < data.nodes.size(); ++i) {
      PointSample &p = data.nodes[i];
      if (p.x != i % nx || p.y != (i / nx) % ny || p.z != i / (nx * ny))
        throw std::runtime_error("Unexpected Cartesian source ordering.");
      p.x = origin[0] + p.x * spacing[0];
      p.y = origin[1] + p.y * spacing[1];
      p.z = origin[2] + p.z * spacing[2];
    }
  }
}
void write_value(std::ostream &out, double v) {
  if (without_sgems::finite(v))
    out << v;
  else
    out << "nan";
}
void run_config(bool templ, const std::string &directory, const std::string &filename) {
  try {
    const Config c(directory, filename);
    const Source source(c);
    std::string out = c.path(c.get("runtime", "output_dir", "output"));
    if (_mkdir(out.c_str()) != 0 && errno != EEXIST)
      throw std::runtime_error("Cannot create result directory.");
    if (templ)
      calculate_templates(c, source, out);
    else
      calculate_windows(c, source, out);
  } catch (const std::bad_alloc &) {
    throw std::runtime_error(
        "Insufficient memory inside 32-bit SGeMS. Reduce the input or template workload; the "
        "standalone 64-bit calculation remains available for larger jobs.");
  }
}

std::vector<double> point_weights(const Config &c, const PointSetData &data, std::ostream &meta) {
  using namespace without_sgems;
  std::string mode = c.get("weighted_moment", "weight_mode", "uniform");
  if (mode == "external") {
    std::vector<double> w(data.nodes.size(), 0);
    size_t count = 0;
    for (size_t i = 0; i < w.size(); ++i) {
      const PointSample& p = data.nodes[i];
      if (!p.informed) continue;
      if (!finite(p.external_weight) || p.external_weight <= 0) {
        std::ostringstream message;
        message << "External weight must be finite and strictly positive at source record " << p.index;
        throw std::runtime_error(message.str());
      }
      w[i] = p.external_weight;
      ++count;
    }
    meta << "weight_mode=external\nweight_source=external_observation_property\nweight_property="
         << c.get("weighted_moment", "source_weight_property", c.get("weighted_moment", "weight_property"))
         << "\nweight_source_object=" << c.get("weighted_moment", "source_weight_object")
         << "\nweight_source_file=" << c.get("general", "source_grid")
         << "\nweight_record_mapping=source_node_id\nweight_valid_count=" << count << '\n';
    return w;
  }
  int dims = integer(number(c.get("declustering", "dimensions", data.dimensions == 2 ? "2" : "3")),
                     2, data.dimensions);
  meta << "weight_mode=" << mode << "\ndeclustering_dimensions=" << dims << '\n';
  if (mode == "uniform") {
    std::vector<double> w(data.nodes.size(), 0);
    for (size_t i = 0; i < w.size(); ++i)
      if (data.nodes[i].informed)
        w[i] = 1;
    return w;
  }
  if (mode == "cell_declustering") {
    CellDeclusteringOptions options;
    options.dimensions = dims;
    std::vector<double> v = c.list("declustering", "cell_size");
    if (v.size() != (size_t)dims && !(dims == 2 && v.size() == 3))
      throw std::runtime_error("Cell size must specify every active axis.");
    for (size_t k = 0; k < v.size(); ++k) {
      if (v[k] <= 0)
        throw std::runtime_error("Cell sizes must be positive.");
      options.cell_size[k] = v[k];
    }
    if (!c.get("declustering", "origin").empty()) {
      v = c.list("declustering", "origin");
      if (v.size() != (size_t)dims && !(dims == 2 && v.size() == 3))
        throw std::runtime_error("Invalid cell origin dimensions.");
      if (dims == 2 && v.size() == 3 && v[2] != 0)
        throw std::runtime_error("2D declustering origin requires Z=0.");
      for (size_t k = 0; k < v.size(); ++k)
        options.origin[k] = v[k];
      options.origin_at_data_min = false;
    }
    options.shifts_per_axis = integer(number(c.get("declustering", "shifts_per_axis", "3")), 1);
    std::string offsets = c.get("declustering", "origin_offsets_file");
    if (!offsets.empty()) {
      if (offsets.size() < 2 || offsets[1] != ':')
        throw std::runtime_error("Origin offsets file requires an absolute path.");
      std::ifstream in(offsets.c_str());
      if (!in)
        throw std::runtime_error("Cannot open origin offsets file.");
      std::string line;
      while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#')
          continue;
        v = numbers(line);
        if (v.size() != (size_t)dims)
          throw std::runtime_error("Origin offset row dimensions do not match.");
        Vec3 row;
        for (int k = 0; k < dims; ++k)
          row[k] = v[k];
        options.origin_offsets.push_back(row);
      }
      if (options.origin_offsets.empty())
        throw std::runtime_error("Origin offsets file is empty.");
    }
    CellDeclusteringResult r = calculate_cell_declustering(data, options);
    meta << "origin_count=" << r.origin_count << "\nweighted_mean=" << r.weighted_mean << '\n';
    return r.weights;
  }
  DeclusteringDomain domain;
  domain.dimensions = dims;
  std::vector<double> lo = c.list("declustering", "domain_min"),
                      hi = c.list("declustering", "domain_max");
  if (lo.size() != (size_t)dims || hi.size() != lo.size())
    throw std::runtime_error("Domain bounds must specify every active axis.");
  for (int k = 0; k < dims; ++k) {
    domain.lower[k] = lo[k];
    domain.upper[k] = hi[k];
  }
  SpatialDeclusteringResult r;
  if (mode == "polygonal_declustering")
    r = calculate_polygonal_declustering(data, domain);
  else if (mode == "kriging_declustering") {
    KrigingDeclusteringOptions options;
    options.domain = domain;
    options.covariance_model = lower(c.get("declustering", "kriging_model"));
    options.nugget_fraction = number(c.get("declustering", "kriging_nugget_fraction"));
    std::vector<double> ranges = c.list("declustering", "kriging_ranges"),
                        cells = c.list("declustering", "kriging_integration_cells");
    if (ranges.size() != (size_t)dims || cells.size() != ranges.size())
      throw std::runtime_error(
          "Kriging ranges and integration counts must specify every active axis.");
    for (int k = 0; k < dims; ++k) {
      options.ranges[k] = ranges[k];
      options.integration_cells[k] = integer(cells[k], 1);
    }
    r = calculate_kriging_declustering(data, options);
  } else
    throw std::runtime_error("Unsupported template weight method.");
  meta << "unique_sites=" << r.unique_site_count << "\nweighted_mean=" << r.weighted_mean
       << "\npartition_error=" << r.partition_relative_error
       << "\nkriging_residual=" << r.kriging_max_residual << '\n';
  return r.weights;
}
} // namespace hosc_native

namespace without_sgems {
namespace {
double coord(const PointSample &p, int a) { return a == 0 ? p.x : a == 1 ? p.y : p.z; }
struct AxisLess {
  const PointSetData &data;
  int axis;
  AxisLess(const PointSetData &d, int a) : data(d), axis(a) {}
  bool operator()(size_t a, size_t b) const {
    double x = coord(data.nodes[a], axis), y = coord(data.nodes[b], axis);
    return x < y || (x == y && a < b);
  }
};
} // namespace
PointSpatialIndex::PointSpatialIndex(const PointSetData &d) : data_(d), root_(-1) {
  std::vector<size_t> ids;
  for (size_t i = 0; i < d.nodes.size(); ++i)
    if (d.nodes[i].informed)
      ids.push_back(i);
  branches_.reserve(ids.size());
  root_ = build(ids, 0, ids.size(), 0);
}
int PointSpatialIndex::build(std::vector<size_t> &ids, size_t lo, size_t hi, int depth) {
  if (lo >= hi)
    return -1;
  int axis = depth % 3;
  size_t mid = lo + (hi - lo) / 2;
  std::nth_element(ids.begin() + lo, ids.begin() + mid, ids.begin() + hi, AxisLess(data_, axis));
  int index = (int)branches_.size();
  branches_.push_back(Branch(ids[mid], axis));
  int l = build(ids, lo, mid, depth + 1), r = build(ids, mid + 1, hi, depth + 1);
  branches_[index].left = l;
  branches_[index].right = r;
  return index;
}
void PointSpatialIndex::query(int index, const Vec3 &lo, const Vec3 &hi,
                              std::vector<size_t> &out) const {
  if (index < 0)
    return;
  const Branch &b = branches_[index];
  const PointSample &p = data_.nodes[b.point];
  if (p.x >= lo[0] && p.x <= hi[0] && p.y >= lo[1] && p.y <= hi[1] && p.z >= lo[2] && p.z <= hi[2])
    out.push_back(b.point);
  double v = coord(p, b.axis);
  if (lo[b.axis] <= v)
    query(b.left, lo, hi, out);
  if (hi[b.axis] >= v)
    query(b.right, lo, hi, out);
}
std::vector<size_t> PointSpatialIndex::query_box(const Vec3 &lo, const Vec3 &hi) const {
  std::vector<size_t> ids;
  query(root_, lo, hi, ids);
  std::sort(ids.begin(), ids.end());
  return ids;
}
std::vector<size_t> PointSpatialIndex::query_radius(const Vec3 &c, double r) const {
  std::vector<size_t> candidates = query_box(Vec3(c[0] - r, c[1] - r, c[2] - r),
                                             Vec3(c[0] + r, c[1] + r, c[2] + r)),
                      out;
  for (size_t i = 0; i < candidates.size(); ++i) {
    const PointSample &p = data_.nodes[candidates[i]];
    if (hypot3(p.x - c[0], p.y - c[1], p.z - c[2]) <= r)
      out.push_back(candidates[i]);
  }
  return out;
}
} // namespace without_sgems
