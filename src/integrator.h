#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <complex>
#include <vector>

#include "math_utils.h"

class PredictorCorrector {
 public:
  PredictorCorrector(const int, const int, const double, const double);

  Eigen::ArrayXXd predictor_coefs;
  Eigen::ArrayXXd corrector_coefs;
  double future_coef;

  size_t width() const { return n_time; }
 private:
  size_t n_time;
};

#endif
