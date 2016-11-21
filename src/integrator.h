#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <complex>
#include <vector>

#include "math_utils.h"

namespace PredictorCorrector {
  class Weights;
}

class PredictorCorrector::Weights {
 public:
  Weights(const size_t, const size_t, const double, const double);

  Eigen::ArrayXXd predictor_coefs;
  Eigen::ArrayXXd corrector_coefs;
  double future_coef;

  size_t width() const { return n_time_; }
 private:
  size_t n_time_;
};

#endif
