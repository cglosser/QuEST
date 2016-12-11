#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <Eigen/Dense>
#include <algorithm>
#include <boost/multi_array.hpp>
#include <complex>
#include <vector>

#include "math_utils.h"
#include "quantum_dot.h"

namespace PredictorCorrector {
  class Weights;
  class Integrator;
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
             const double);

  void step();

 private:
  typedef boost::multi_array<matrix_elements, 3,
                             Eigen::aligned_allocator<matrix_elements>>
      HistoryArray;

  int now;
  int num_solutions, num_steps;
  double dt;
  Weights weights;
  HistoryArray history;

  void predictor();
  void evaluator();
  void corrector();
};

#endif
