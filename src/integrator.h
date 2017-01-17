#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <Eigen/Dense>
#include <algorithm>
#include <boost/multi_array.hpp>
#include <complex>
#include <vector>

#include "history.h"
#include "math_utils.h"

namespace PredictorCorrector {
  class Weights;
  class Integrator;

  typedef std::function<History::soltype(const History::soltype &,
                                         const double)>
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
  Integrator(const int, const int, const double, const int, const int,
             const double, const std::shared_ptr<History::HistoryArray>,
             const std::vector<rhs_func> &);
  void solve(const int);

 private:
  int num_solutions, num_steps;
  double dt;
  Weights weights;
  std::shared_ptr<History::HistoryArray> history;
  std::vector<rhs_func> rhs_funcs;

  void predictor(const int);
  void evaluator(const int);
  void corrector(const int);
};

#endif
