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
  Weights(const int, const int, const double, const double);

  Eigen::ArrayXXd ps, cs;
  double future_coef;

  int width() const { return n_time_; }
 private:
  int n_time_;
};

class PredictorCorrector::Integrator {
 public:
  Integrator(std::vector<QuantumDot>, const int, const double, const int,
             const int, const double, const double);

  void step();

  // private:

  typedef boost::multi_array<matrix_elements, 3,
                             Eigen::aligned_allocator<matrix_elements>>
      HistoryArray;

  int now;
  int num_steps;
  double dt;
  Weights weights;
  HistoryArray history;

  std::vector<QuantumDot> dots;
  std::vector<double> rabi_freqs;

  void predictor();
  void evaluator();
  void corrector();

  void set_rabi_freqs(const double);
};

#endif
