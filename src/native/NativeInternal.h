#ifndef HOSC_NATIVE_INTERNAL_H
#define HOSC_NATIVE_INTERNAL_H
#include "NativeTypes.h"
#include <climits>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>

namespace hosc_native {
using without_sgems::PointSample;
using without_sgems::PointSetData;
using without_sgems::Vec3;
std::string trim(const std::string &);
std::string lower(std::string);
double number(const std::string &);
std::vector<double> numbers(const std::string &);
int integer(double, int minimum = 0, int maximum = INT_MAX);
Vec3 vector3(const std::vector<double> &);
struct Config {
  std::map<std::string, std::string> items;
  std::string directory;
  Config(const std::string &, const std::string &);
  std::string get(const std::string &, const std::string &, const std::string &fallback = "") const;
  std::vector<double> list(const std::string &, const std::string &,
                           const std::string &fallback = "") const;
  bool flag(const std::string &, const std::string &, bool fallback = false) const;
  std::string path(const std::string &) const;
};
struct Source {
  PointSetData data;
  bool points;
  int nx, ny, nz;
  Vec3 origin, spacing;
  Source(const Config &);
};
typedef std::vector<std::vector<double> > Values;
typedef std::vector<std::vector<size_t> > Ids;
double evaluate(Values &, const std::vector<double> &, const std::vector<int> &, bool moment,
                bool bell);
std::vector<double> point_weights(const Config &, const PointSetData &, std::ostream &);
void calculate_templates(const Config &, const Source &, const std::string &);
void calculate_windows(const Config &, const Source &, const std::string &);
void write_value(std::ostream &, double);
} // namespace hosc_native
#endif
