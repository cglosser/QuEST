#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <Eigen/Dense>
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
  Weights(const size_t, const size_t, const double, const double);

  Eigen::ArrayXXd ps, cs;
  double future_coef;

  size_t width() const { return n_time_; }
 private:
  size_t n_time_;
};

class PredictorCorrector::Integrator {
 public:
  Integrator(std::vector<QuantumDot>, const size_t, const double, const size_t,
             const size_t, const double, const double);

  // private:

  typedef boost::multi_array<matrix_elements, 3,
                             Eigen::aligned_allocator<matrix_elements>>
      HistoryArray;

  int now;
  size_t num_steps;
  double dt;
  Weights weights;
  HistoryArray history;
  std::vector<QuantumDot> dots;
};

#endif
