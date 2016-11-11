#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <complex>
#include <cmath>
#include <iostream>
#include <vector>

#include "math_utils.h"

class PredictorCorrector {
 public:
  PredictorCorrector(const int, const int, const double, const double);

  size_t n_lambda, n_time;
  double rho, tolerance, timestep, future_time, step_factor;
  Eigen::VectorXcd lambdas;
  Eigen::VectorXd times;

 //private:

  Eigen::MatrixXcd predictor_matrix() const;
  Eigen::MatrixXcd corrector_matrix() const;

};

#endif
