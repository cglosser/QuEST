#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <Eigen/Dense>
#include <algorithm>
#include <boost/multi_array.hpp>
#include <complex>
#include <vector>

#include "common.h"
#include "interaction_table.h"
#include "math_utils.h"

namespace PredictorCorrector {
  class Weights;
  class Integrator;

  typedef std::function<soltype(const soltype &, const double)> rhs_func;
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
             const double, const std::vector<rhs_func> &, InteractionTable &);
  void step();

 private:
  int now;
  int num_solutions, num_steps;
  double dt;
  Weights weights;
  HistoryArray history;
  std::vector<rhs_func> rhs_funcs;
  InteractionTable interaction_table;

  void predictor();
  void evaluator();
  void corrector();
};

#endif
