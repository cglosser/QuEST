#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <Eigen/Dense>
#include <algorithm>
#include <boost/multi_array.hpp>
#include <complex>
#include <string>
#include <vector>

#include "history.h"
#include "interactions/pulse_interaction.h"
#include "math_utils.h"

namespace PredictorCorrector {
  class Weights;
  class Integrator;

  typedef std::function<History::soltype(const History::soltype &, const cmplx)>
      rhs_func;
}

class PredictorCorrector::Weights {
 public:
  Weights(const int, const int, const double);

  Eigen::ArrayXXd ps, cs;
  double future_coef;

  int width() const { return n_time; }
 private:
  int n_time;
};

class PredictorCorrector::Integrator {
 public:
  Integrator(const double, const int, const int, const double,
             const std::shared_ptr<History::HistoryArray> &,
             const std::vector<rhs_func> &,
             std::vector<std::shared_ptr<Interaction>>);
  void solve() const;

 private:
  int num_solutions, time_idx_ubound;
  double dt;
  Weights weights;
  std::shared_ptr<History::HistoryArray> history;
  std::vector<rhs_func> rhs_funcs;
  std::vector<std::shared_ptr<Interaction>> interactions;

  void solve_step(const int) const;
  void predictor(const int) const;
  void evaluator(const int) const;
  void corrector(const int) const;

  bool all_finite(const int) const;
  void throw_if_unbounded_solution(const int) const;
  void log_percentage_complete(const int) const;
};

#endif
