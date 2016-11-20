#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <complex>
#include <vector>

#include "math_utils.h"

class PredictorCorrector {
 public:
  PredictorCorrector(const int, const int, const double, const double);

  double future_coef;
  Eigen::ArrayXXd predictor_coefs;
  Eigen::ArrayXXd corrector_coefs;

  size_t num_steps() const { return n_time; }

 private:
  size_t n_lambda, n_time;
  double step_factor, rho, tolerance, timestep, future_time;
  Eigen::VectorXcd lambdas;
  Eigen::VectorXd times;

  Eigen::MatrixXcd predictor_matrix() const;
  Eigen::MatrixXcd corrector_matrix() const;
  Eigen::VectorXcd rhs_vector() const;
  Eigen::VectorXd compute_coefficients(const Eigen::MatrixXcd &) const;
};

#endif
