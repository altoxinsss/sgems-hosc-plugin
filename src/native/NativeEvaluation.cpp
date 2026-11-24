#include "NativeInternal.h"
#include "WeightedMomentCore.h"
#include "vendor/FunctionCall.h"
#include <new>

namespace hosc_native {
namespace wmc = weighted_moment_core;
namespace {
class Workspace {
  int variables_;

public:
  Workspace(int n) : variables_(n) { ClearOneEvaluationWorkspaceForDLL(n); }
  ~Workspace() {
    SetLowOrderCumulantMomentCacheMode(0);
    ResetFactorSwitchesForRun();
    ClearOneEvaluationWorkspaceForDLL(variables_);
    G_Template = NULL;
    G_TemplateNum = 0;
  }
};
void store_moment(size_t slot, const std::vector<int> &powers,
                  const std::vector<wmc::WeightedSample> &samples,
                  const wmc::WeightedMomentOptions &options) {
  int variables = (int)powers.size() - 1;
  G_Moment[slot] =
      wmc::weighted_moment(&samples[0], (int)samples.size(), &powers[1], variables, options);
  if (!without_sgems::finite(G_Moment[slot]))
    throw std::runtime_error("Template moment overflow.");
  G_MomentIndex[slot] = (char *)calloc(variables + 1, sizeof(char));
  if (!G_MomentIndex[slot])
    throw std::bad_alloc();
  int last = variables;
  while (last > 0 && powers[last] == 0)
    --last;
  for (int j = 1; j <= last; ++j)
    G_MomentIndex[slot][j - 1] = (char)('0' + powers[j]);
}
void supplied_moments(const std::vector<wmc::WeightedSample> &samples,
                      const std::vector<int> &index, const wmc::WeightedMomentOptions &options) {
  int variables = (int)index.size() - 1;
  size_t count = 1;
  for (int j = 1; j <= variables; ++j)
    count *= index[j] + 1;
  G_SeriesNum = (int)count;
  G_Moment = (double *)calloc(count, sizeof(double));
  G_MomentIndex = (char **)calloc(count, sizeof(char *));
  if (!G_Moment || !G_MomentIndex)
    throw std::bad_alloc();
  store_moment(0, index, samples, options);
  for (size_t code = 0; code + 1 < count; ++code) {
    size_t remainder = code;
    std::vector<int> powers(variables + 1, 0);
    for (int j = variables; j >= 1; --j) {
      powers[j] = (int)(remainder % (index[j] + 1));
      powers[0] += powers[j];
      remainder /= index[j] + 1;
    }
    store_moment(code + 1, powers, samples, options);
  }
  G_CurrMomentPos = G_SeriesNum;
  G_MomentFull = true;
  G_MomentHashReady = false;
  if (BuildSuppliedMomentHashForDLL() != 0)
    throw std::runtime_error("Cannot build moment lookup index.");
}
} // namespace
double evaluate(Values &values, const std::vector<double> &weights, const std::vector<int> &powers,
                bool moment, bool bell) {
  if (values.empty())
    return without_sgems::nan();
  if (values.size() > (size_t)INT_MAX)
    throw std::runtime_error("Sample count exceeds evaluator integer range.");
  int variables = (int)powers.size(), order = 0;
  for (int j = 0; j < variables; ++j)
    order += powers[j];
  if (order < 1)
    throw std::runtime_error("Order must be positive.");
  bool weighted = !weights.empty();
  if (weighted && weights.size() != values.size())
    throw std::runtime_error("Tuple weight count mismatch.");
  wmc::WeightedMomentOptions options;
  options.mode = weighted ? wmc::WeightMode_External : wmc::WeightMode_Uniform;
  std::vector<wmc::WeightedSample> samples(values.size());
  for (size_t i = 0; i < values.size(); ++i) {
    samples[i].values = &values[i][0];
    samples[i].value_count = variables;
    samples[i].weight = weighted ? weights[i] : 1;
  }
  if (moment || order == 1)
    return wmc::weighted_moment(&samples[0], (int)samples.size(), &powers[0], variables, options);
  for (int j = 0; j < variables; ++j) {
    long double sum = 0, total = 0;
    for (size_t i = 0; i < values.size(); ++i) {
      long double w = weighted ? weights[i] : 1;
      sum += w * values[i][j];
      total += w;
    }
    long double mean = sum / total;
    for (size_t i = 0; i < values.size(); ++i)
      values[i][j] = (double)(values[i][j] - mean);
  }
  if (bell && variables == 1) {
    std::vector<double> v(values.size());
    for (size_t i = 0; i < v.size(); ++i)
      v[i] = values[i][0];
    return wmc::weighted_univariate_cumulant(&v[0], weighted ? &weights[0] : NULL, (int)v.size(),
                                             order);
  }
  if (bell && weighted)
    return wmc::weighted_cumulant_from_multiplicities(&samples[0], (int)samples.size(), &powers[0],
                                                      variables, options);
  std::vector<int> index(variables + 1);
  index[0] = order;
  for (int j = 0; j < variables; ++j)
    index[j + 1] = powers[j];
  Values storage(values.size(), std::vector<double>(variables + 1, 0));
  std::vector<Template> views(values.size());
  for (size_t i = 0; i < views.size(); ++i) {
    std::copy(values[i].begin(), values[i].end(), storage[i].begin() + 1);
    views[i].Value = &storage[i][0];
    views[i].Lag = NULL;
    views[i].Mean = 0;
    views[i].index = &index[0];
  }
  Workspace workspace(variables);
  G_Template = &views[0];
  G_TemplateNum = (int)views.size();
  double result = myEmpty;
  if (bell)
    result = BellCumulantFromMultiplicityIndex(&index[0], variables);
  else {
    SetFactorSwitchesForRun(1, 1, 1, 1);
    SetLowOrderCumulantMomentCacheMode(1);
    if (weighted)
      supplied_moments(samples, index, options);
    if ((weighted || PrecomputeMomentsForIndex(&index[0], variables) != myEmpty) &&
        InitAVLWorkspaceForDLL(variables) == 0)
      result = BetterRecCumulant3(&index[0], variables);
  }
  if (result == myEmpty || !without_sgems::finite(result))
    throw std::runtime_error("Bell/AVL evaluation failed or overflowed.");
  return result;
}
} // namespace hosc_native
