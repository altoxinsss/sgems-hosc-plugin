#include "NativeInternal.h"
#include "vendor/FunctionCall.h"
#include <new>

namespace hosc_native {
namespace {
const double pi = 3.14159265358979323846;
struct Direction {
  int nod, count;
  double lag, ltol, azm, atol, bandwh, dip, dtol, bandwd;
  Vec3 offset;
  Direction()
      : nod(1), count(1), lag(1), ltol(0), azm(0), atol(0), bandwh(0), dip(0), dtol(0), bandwd(0) {}
};
std::vector<Direction> read_directions(const Config &c, bool points, bool exact) {
  const char *names[] = {"nod",  "nlag",   "lag", "ltol", "azm",
                         "atol", "bandwh", "dip", "dtol", "bandwd"};
  std::vector<double> lists[10];
  for (int k = 0; k < 10; ++k)
    lists[k] = c.list("template_shape", names[k]);
  size_t n = lists[0].size();
  if (n < 1 || n > 8)
    throw std::runtime_error("Template needs one to eight directions.");
  for (int k = 1; k < 10; ++k)
    if (lists[k].size() != n)
      throw std::runtime_error("Directional lists must have matching lengths.");
  std::vector<Direction> out(n);
  int order = 1;
  for (size_t j = 0; j < n; ++j) {
    Direction &d = out[j];
    d.nod = integer(lists[0][j], 1, 8);
    d.count = integer(lists[1][j], 1, 100000);
    d.lag = lists[2][j];
    d.ltol = lists[3][j];
    d.azm = lists[4][j];
    d.atol = lists[5][j];
    d.bandwh = lists[6][j];
    d.dip = lists[7][j];
    d.dtol = lists[8][j];
    d.bandwd = lists[9][j];
    order += d.nod;
    if (d.lag == 0 || d.ltol < 0 || d.atol < 0 || d.atol > 180 || d.bandwh < 0 ||
        fabs(d.dip) > 90 || d.dtol < 0 || d.dtol > 180 || d.bandwd < 0)
      throw std::runtime_error("Invalid lag, tolerance, bandwidth or dip.");
    if (!points) {
      for (int k = 0; k < 10; ++k)
        integer(lists[k][j], INT_MIN);
      integer(d.lag, 1);
      if (!exact && d.atol >= 90)
        throw std::runtime_error("legacy_search requires ATOL<90 degrees.");
      if ((double)(d.count - 1) * d.lag > INT_MAX)
        throw std::runtime_error("Cartesian lag distance exceeds integer range.");
    }
  }
  if (exact) {
    std::string text = c.get("template_shape", "integer_direction_vectors");
    size_t start = 0, j = 0;
    while (start <= text.size()) {
      size_t end = text.find(';', start);
      std::vector<double> v =
          numbers(text.substr(start, end == std::string::npos ? end : end - start));
      if (j >= n || v.size() < 2 || v.size() > 3)
        throw std::runtime_error("Invalid integer direction vectors.");
      for (size_t k = 0; k < v.size(); ++k)
        out[j].offset[k] = integer(v[k], INT_MIN);
      if (out[j].offset[0] == 0 && out[j].offset[1] == 0 && out[j].offset[2] == 0)
        throw std::runtime_error("Integer direction vector cannot be zero.");
      for (int k = 0; k < 3; ++k)
        if (fabs(out[j].offset[k]) * (out[j].count - 1) * out[j].lag > INT_MAX)
          throw std::runtime_error("Template offset exceeds integer range.");
      ++j;
      if (end == std::string::npos)
        break;
      start = end + 1;
    }
    if (j != n)
      throw std::runtime_error("One integer vector is required per direction.");
  }
  return out;
}
Vec3 displacement(const Direction &d, int p) {
  double az = d.azm * pi / 180, dip = d.dip * pi / 180, lag = p * d.lag;
  return Vec3(lag * cos(dip) * cos(az), lag * cos(dip) * sin(az), lag * sin(dip));
}
bool in_direction(const PointSample &a, const PointSample &p, const Direction &d, int position,
                  int dimensions) {
  double dx = p.x - a.x, dy = p.y - a.y, dz = p.z - a.z, dist = without_sgems::hypot3(dx, dy, dz),
         signedlag = position * d.lag, lag = fabs(signedlag),
         eps = 64 * DBL_EPSILON * std::max(1.0, std::max(lag, dist));
  if (dist <= eps || fabs(dist - lag) > d.ltol + eps)
    return false;
  double az = d.azm + (signedlag < 0 ? 180 : 0), dip = signedlag < 0 ? -d.dip : d.dip,
         horiz = _hypot(dx, dy);
  double da = fmod(atan2(dy, dx) * 180 / pi - az, 360.0);
  if (da > 180)
    da -= 360;
  if (da < -180)
    da += 360;
  if (horiz > eps && fabs(cos(dip * pi / 180)) > 1e-12 && fabs(da) > d.atol + 1e-10)
    return false;
  if (dimensions == 3 && fabs(atan2(dz, horiz) * 180 / pi - dip) > d.dtol + 1e-10)
    return false;
  double across = -dx * sin(az * pi / 180) + dy * cos(az * pi / 180),
         vertical = -sin(dip * pi / 180) * (dx * cos(az * pi / 180) + dy * sin(az * pi / 180)) +
                    cos(dip * pi / 180) * dz;
  return (d.bandwh == 0 || fabs(across) <= d.bandwh + eps) &&
         (d.bandwd == 0 || fabs(vertical) <= d.bandwd + eps);
}
struct Samples {
  Values values;
  Ids ids;
  size_t anchors;
  Samples() : anchors(0) {}
};
struct Enumerate {
  const PointSetData &data;
  const Ids &candidates;
  const std::vector<int> &same;
  std::vector<size_t> selected;
  Samples &out;
  size_t max_tuples, max_visits;
  size_t &visits;
  Enumerate(const PointSetData &d, const Ids &c, const std::vector<int> &s, size_t anchor,
            Samples &o, size_t mt, size_t mv, size_t &v)
      : data(d), candidates(c), same(s), selected(c.size(), anchor), out(o), max_tuples(mt),
        max_visits(mv), visits(v) {}
  void visit(size_t j) {
    if (max_visits) {
      if (visits >= max_visits)
        throw std::runtime_error("Template candidate enumeration exceeds max_candidate_visits; no "
                                 "tuples were truncated.");
      ++visits;
    }
    if (j == selected.size()) {
      if (max_tuples && out.values.size() >= max_tuples)
        throw std::runtime_error("Template exceeds max_template_tuples; no tuples were truncated.");
      if (out.values.size() >= (size_t)INT_MAX)
        throw std::runtime_error("Template count exceeds evaluator integer range.");
      std::vector<double> row(selected.size());
      for (size_t k = 0; k < row.size(); ++k)
        row[k] = data.nodes[selected[k]].value;
      out.values.push_back(row);
      out.ids.push_back(selected);
      return;
    }
    if (same[j] >= 0) {
      size_t id = selected[same[j]];
      if (std::binary_search(candidates[j].begin(), candidates[j].end(), id)) {
        selected[j] = id;
        visit(j + 1);
      }
      return;
    }
    for (size_t i = 0; i < candidates[j].size(); ++i) {
      size_t id = candidates[j][i];
      if (std::find(selected.begin(), selected.begin() + j, id) != selected.begin() + j)
        continue;
      selected[j] = id;
      visit(j + 1);
    }
  }
};
size_t count_limit(const Config &c, const char *key) {
  std::string s = c.get("template_shape", key, "0");
  if (s.empty() || s.find_first_not_of("0123456789") != std::string::npos)
    throw std::runtime_error("Count limit must be an unsigned integer.");
  unsigned __int64 v = 0;
  for (size_t i = 0; i < s.size(); ++i) {
    v = v * 10 + s[i] - '0';
    if (v > UINT_MAX)
      throw std::runtime_error("Count limit exceeds 32-bit addressable range.");
  }
  return (size_t)v;
}
Samples collect_points(const Source &source, const without_sgems::PointSpatialIndex &spatial,
                       const std::vector<Direction> &dirs, const std::vector<int> &pos, size_t mt,
                       size_t mv) {
  const PointSetData &data = source.data;
  size_t vars = dirs.size() + 1;
  std::vector<Vec3> targets(vars);
  std::vector<int> same(vars, -1);
  Samples out;
  for (size_t j = 1; j < vars; ++j) {
    targets[j] = displacement(dirs[j - 1], pos[j - 1]);
    double scale = 1;
    for (int axis = 0; axis < 3; ++axis)
      scale = std::max(scale, fabs(targets[j][axis]));
    for (size_t k = 0; k < j; ++k)
      if (without_sgems::hypot3(targets[j][0] - targets[k][0], targets[j][1] - targets[k][1],
                                targets[j][2] - targets[k][2]) <= 64 * DBL_EPSILON * scale) {
        same[j] = (int)k;
        break;
      }
  }
  size_t visits = 0;
  for (size_t anchor = 0; anchor < data.nodes.size(); ++anchor) {
    const PointSample &a = data.nodes[anchor];
    if (!a.informed)
      continue;
    Ids candidates(vars);
    candidates[0].push_back(anchor);
    bool possible = true;
    for (size_t j = 1; j < vars; ++j) {
      if (pos[j - 1] == 0)
        candidates[j].push_back(anchor);
      else {
        const Direction &d = dirs[j - 1];
        double radius = fabs(pos[j - 1] * d.lag) + d.ltol;
        radius += 64 * DBL_EPSILON * std::max(1.0, radius);
        std::vector<size_t> nearby = spatial.query_radius(Vec3(a.x, a.y, a.z), radius);
        for (size_t k = 0; k < nearby.size(); ++k)
          if (nearby[k] != anchor &&
              in_direction(a, data.nodes[nearby[k]], d, pos[j - 1], data.dimensions))
            candidates[j].push_back(nearby[k]);
      }
      if (candidates[j].empty()) {
        possible = false;
        break;
      }
    }
    if (!possible)
      continue;
    size_t before = out.values.size();
    Enumerate task(data, candidates, same, anchor, out, mt, mv, visits);
    task.visit(1);
    if (out.values.size() > before)
      ++out.anchors;
  }
  return out;
}
Samples collect_exact(const Source &s, const std::vector<Direction> &dirs,
                      const std::vector<int> &pos) {
  Samples out;
  for (size_t a = 0; a < s.data.nodes.size(); ++a) {
    if (!s.data.nodes[a].informed)
      continue;
    std::vector<size_t> ids(1, a);
    std::vector<double> row(1, s.data.nodes[a].value);
    bool complete = true;
    for (size_t j = 0; j < dirs.size(); ++j) {
      double distance = pos[j] * dirs[j].lag;
      __int64 x = (int)(a % s.nx) + (__int64)(dirs[j].offset[0] * distance),
              y = (int)((a / s.nx) % s.ny) + (__int64)(dirs[j].offset[1] * distance),
              z = (int)(a / (s.nx * s.ny)) + (__int64)(dirs[j].offset[2] * distance);
      if (x < 0 || y < 0 || z < 0 || x >= s.nx || y >= s.ny || z >= s.nz) {
        complete = false;
        break;
      }
      size_t id = ((size_t)z * s.ny + (size_t)y) * s.nx + (size_t)x;
      if (!s.data.nodes[id].informed) {
        complete = false;
        break;
      }
      ids.push_back(id);
      row.push_back(s.data.nodes[id].value);
    }
    if (complete) {
      out.values.push_back(row);
      out.ids.push_back(ids);
      ++out.anchors;
    }
  }
  return out;
}
struct Trace {
  const Source *source;
  const std::vector<int> *order;
  Samples *samples;
};
void trace(int vars, const int *x, const int *y, const int *z, const double *, void *context) {
  Trace &t = *(Trace *)context;
  std::vector<size_t> ids(vars);
  std::vector<double> row(vars);
  for (int j = 0; j < vars; ++j) {
    size_t slot = j == 0 ? 0 : (*t.order)[j - 1] + 1;
    size_t id = ((size_t)z[j] * t.source->ny + y[j]) * t.source->nx + x[j];
    ids[slot] = id;
    row[slot] = t.source->data.nodes[id].value;
  }
  t.samples->ids.push_back(ids);
  t.samples->values.push_back(row);
}
class LegacySampler {
  const Source &s_;
  Keys *key_;
  std::vector<int> order_;
  std::vector<unsigned char> informed_;
  bool initialized_;
  void release() {
    SetTemplateRawTrace(NULL, NULL, NULL, 0, 0);
    if (initialized_) {
      if (s_.nz == 1) {
        Free2DVec(&ValueMap_2D, MaxX);
        ValueMap_2D = NULL;
      } else {
        Free3DVec(&ValueMap_3D, MaxX, MaxY);
        ValueMap_3D = NULL;
      }
      Destruction();
      initialized_ = false;
    }
    if (key_) {
      int *fields[] = {key_->NOD,  key_->NLAG,   key_->LAG, key_->LTOL, key_->AZM,
                       key_->ATOL, key_->BANDWH, key_->DIP, key_->DTOL, key_->BANDWD};
      for (int k = 0; k < 10; ++k)
        free(fields[k]);
      free(key_);
      key_ = NULL;
      pt_Keys = NULL;
    }
  }

public:
  LegacySampler(const Source &s, const std::vector<Direction> &dirs, const std::vector<int> &order)
      : s_(s), key_(NULL), order_(order), initialized_(false) {
    try {
      key_ = (Keys *)calloc(1, sizeof(Keys));
      if (!key_)
        throw std::bad_alloc();
      key_->DIM = s.nz == 1 ? 2 : 3;
      key_->REG_IRREG = 1;
      key_->NX = s.nx;
      key_->NY = s.ny;
      key_->NZ = s.nz;
      key_->DIMX = s.nx - 1;
      key_->DIMY = s.ny - 1;
      key_->DIMZ = s.nz - 1;
      key_->ORDER = (int)dirs.size() + 1;
      key_->NDIR = (int)dirs.size();
      int **fields[] = {&key_->NOD,  &key_->NLAG,   &key_->LAG, &key_->LTOL, &key_->AZM,
                        &key_->ATOL, &key_->BANDWH, &key_->DIP, &key_->DTOL, &key_->BANDWD};
      for (int k = 0; k < 10; ++k) {
        *fields[k] = (int *)calloc(dirs.size(), sizeof(int));
        if (!*fields[k])
          throw std::bad_alloc();
      }
      for (size_t j = 0; j < dirs.size(); ++j) {
        const Direction &d = dirs[order[j]];
        double v[] = {d.nod,  d.count,  d.lag, d.ltol, d.azm,
                      d.atol, d.bandwh, d.dip, d.dtol, d.bandwd};
        for (int k = 0; k < 10; ++k)
          (*fields[k])[j] = (int)v[k];
      }
      size_t n = s.data.nodes.size();
      informed_.resize(n);
      double **grid = new double *[n + 1];
      grid[0] = new double[2];
      grid[0][0] = (double)n;
      grid[0][1] = 4;
      for (size_t i = 0; i < n; ++i) {
        grid[i + 1] = new double[4];
        grid[i + 1][0] = (double)(i % s.nx);
        grid[i + 1][1] = (double)((i / s.nx) % s.ny);
        grid[i + 1][2] = (double)(i / (s.nx * s.ny));
        grid[i + 1][3] = s.data.nodes[i].informed ? s.data.nodes[i].value : 0;
        informed_[i] = s.data.nodes[i].informed ? 1 : 0;
      }
      if (InitialStruct4(key_, grid) != 0)
        throw std::runtime_error("Cannot initialize legacy Cartesian sampler.");
      initialized_ = true;
    } catch (...) {
      release();
      throw;
    }
  }
  ~LegacySampler() { release(); }
  Samples collect(const std::vector<Direction> &dirs, const std::vector<int> &pos) {
    Samples samples;
    Trace t = {&s_, &order_, &samples};
    std::vector<int> distances;
    for (size_t j = 0; j < order_.size(); ++j)
      distances.push_back((int)(pos[order_[j]] * dirs[order_[j]].lag));
    SetTemplateRawTrace(trace, &t, &informed_[0], s_.nx, s_.ny);
    Template *templates = NULL;
    int count = 0;
    int status =
        s_.nz == 1
            ? TemplateConstruct_2D(&distances[0], key_->ORDER, &templates, &count, 0, pt_Grids->ND)
            : TemplateConstruct_3D(&distances[0], key_->ORDER, &templates, &count, 0, pt_Grids->ND);
    SetTemplateRawTrace(NULL, NULL, NULL, 0, 0);
    if (templates) {
      free(templates[0].index);
      templates[0].index = NULL;
      FreeTemplate(&templates);
      G_Template = NULL;
      G_TemplateNum = 0;
    }
    if (status < 0 || count != (int)samples.values.size())
      throw std::runtime_error("Legacy template trace mismatch.");
    samples.anchors = samples.values.size();
    return samples;
  }
};
struct SamplerOwner {
  LegacySampler *p;
  SamplerOwner() : p(NULL) {}
  ~SamplerOwner() { delete p; }
};
} // namespace

void calculate_templates(const Config &c, const Source &source, const std::string &output) {
  bool points = source.points;
  std::string sampling = c.get("template_shape", "template_sampling");
  bool exact = !points && sampling == "exhaustive_integer_offsets";
  if (!points && !exact && sampling != "legacy_search")
    throw std::runtime_error("Unknown Cartesian template sampling.");
  const std::vector<Direction> dirs = read_directions(c, points, exact);
  std::vector<int> powers(1, 1);
  for (size_t j = 0; j < dirs.size(); ++j)
    powers.push_back(dirs[j].nod);
  if (source.data.dimensions == 2)
    for (size_t j = 0; j < dirs.size(); ++j) {
      if (points && dirs[j].dip != 0)
        throw std::runtime_error("Two-dimensional point sets require DIP=0.");
      if (exact && dirs[j].offset[2] != 0)
        throw std::runtime_error("Two-dimensional integer templates require Z offset=0.");
    }
  std::vector<double> disp = c.list("target_cross_sections", "displayed_dirs"),
                      fixed = c.list("target_cross_sections", "fixed_positions");
  if (disp.empty() || disp.size() > 3)
    throw std::runtime_error("Select one to three display directions.");
  std::vector<int> display, hidden, positions(dirs.size(), 0);
  for (size_t i = 0; i < disp.size(); ++i) {
    int d = integer(disp[i], 1, (int)dirs.size()) - 1;
    if (std::find(display.begin(), display.end(), d) != display.end())
      throw std::runtime_error("Duplicate displayed direction.");
    display.push_back(d);
  }
  for (size_t j = 0; j < dirs.size(); ++j)
    if (std::find(display.begin(), display.end(), (int)j) == display.end())
      hidden.push_back((int)j);
  if (hidden.size() != fixed.size())
    throw std::runtime_error("One fixed position is required per hidden direction.");
  for (size_t j = 0; j < hidden.size(); ++j)
    positions[hidden[j]] = integer(fixed[j], 1, dirs[hidden[j]].count) - 1;
  size_t rows = 1;
  for (size_t j = 0; j < display.size(); ++j) {
    size_t n = dirs[display[j]].count;
    if (rows > (size_t)INT_MAX / n)
      throw std::runtime_error("Output count exceeds supported range.");
    rows *= n;
  }
  std::ofstream meta((output + "\\result.meta.txt").c_str());
  meta << std::setprecision(17)
       << "native_vc71_in_process=true\nsource_type=" << (points ? "point_set" : "cartesian")
       << "\ntemplate_sampling=" << sampling << "\n";
  std::vector<double> weights = point_weights(c, source.data, meta);
  bool weighted = c.get("weighted_moment", "weight_mode", "uniform") != "uniform";
  std::ofstream pointout((output + "\\result.point_weights.csv").c_str());
  pointout << "node_id,x,y,z,value,informed,weight\n" << std::setprecision(17);
  for (size_t i = 0; i < source.data.nodes.size(); ++i) {
    const PointSample &p = source.data.nodes[i];
    pointout << p.index << ',' << p.x << ',' << p.y << ',' << p.z << ',';
    write_value(pointout, p.value);
    pointout << ',' << (p.informed ? 1 : 0) << ',' << weights[i] << '\n';
  }
  std::ofstream out((output + "\\result.txt").c_str()),
      support((output + "\\result.support.csv").c_str());
  out << "HOSC result\n4\nx\ny\nz\nresult\n" << std::setprecision(17);
  support << "lag_cell,lag_positions,lag_distances,anchor_count,support_count,value,status,weight_"
             "sum,ESS\n"
          << std::setprecision(17);
  bool moment = lower(c.get("general", "calculate_target")) == "moment",
       bell = lower(c.get("general", "calculate_algorithm")).find("bell") != std::string::npos;
  size_t minimum = integer(number(c.get("template_shape", "minimum_template_count", "1")), 1),
         mt = count_limit(c, "max_template_tuples"), mv = count_limit(c, "max_candidate_visits");
  Vec3 origin = vector3(c.list("general", "target_origin", "0,0,0"));
  without_sgems::PointSpatialIndex spatial(source.data);
  SamplerOwner legacy;
  if (!points && !exact) {
    std::vector<int> reordered = hidden;
    reordered.insert(reordered.end(), display.begin(), display.end());
    legacy.p = new LegacySampler(source, dirs, reordered);
  }
  for (size_t row = 0; row < rows; ++row) {
    size_t rem = row; /* Point reference enumerates original directions; Cartesian reference sorts
                         displayed directions. */
    if (points) {
      for (int j = (int)dirs.size() - 1; j >= 0; --j)
        if (std::find(display.begin(), display.end(), j) != display.end()) {
          positions[j] = (int)(rem % dirs[j].count);
          rem /= dirs[j].count;
        }
    } else {
      for (int j = (int)display.size() - 1; j >= 0; --j) {
        positions[display[j]] = (int)(rem % dirs[display[j]].count);
        rem /= dirs[display[j]].count;
      }
    }
    Samples samples = points  ? collect_points(source, spatial, dirs, positions, mt, mv)
                      : exact ? collect_exact(source, dirs, positions)
                              : legacy.p->collect(dirs, positions);
    std::vector<double> tuple_weights;
    long double sum = 0, squares = 0;
    for (size_t i = 0; i < samples.ids.size(); ++i) {
      std::vector<size_t> ids = samples.ids[i];
      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
      long double logs = 0;
      for (size_t j = 0; j < ids.size(); ++j) {
        double w = weights[ids[j]];
        if (!(w > 0) || !without_sgems::finite(w))
          throw std::runtime_error("Invalid point weight in complete template.");
        logs += std::log((long double)w);
      }
      double w = (double)std::exp(logs / ids.size());
      tuple_weights.push_back(w);
      sum += w;
      squares += (long double)w * w;
    }
    if (!weighted)
      tuple_weights.clear();
    double value = samples.values.size() < minimum
                       ? without_sgems::nan()
                       : evaluate(samples.values, tuple_weights, powers, moment, bell);
    Vec3 coord;
    for (size_t j = 0; j < display.size(); ++j)
      coord[j] = positions[display[j]] * dirs[display[j]].lag;
    for (int k = 0; k < 3; ++k)
      out << coord[k] + (points ? origin[k] : 0) << ' ';
    write_value(out, value);
    out << '\n';
    support << row << ",\"";
    for (size_t j = 0; j < dirs.size(); ++j) {
      if (j)
        support << ';';
      support << positions[j];
    }
    support << "\",\"";
    for (size_t j = 0; j < dirs.size(); ++j) {
      if (j)
        support << ';';
      support << positions[j] * dirs[j].lag;
    }
    support << "\"," << samples.anchors << ',' << samples.values.size() << ',';
    write_value(support, value);
    support << ',' << (without_sgems::finite(value) ? "ok" : "insufficient_support") << ','
            << (double)sum << ',' << (squares > 0 ? (double)(sum * sum / squares) : 0) << '\n';
  }
  if (!out || !support || !meta || !pointout)
    throw std::runtime_error("Cannot write template outputs.");
}
} // namespace hosc_native
