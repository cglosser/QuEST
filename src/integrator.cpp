#include "integrator.h"

std::complex<double> semidisk(const double);

PredictorCorrector::PredictorCorrector(const int nlam, const int ntim,
  const double radius, const double toler)
    : n_lambda(nlam),
      n_time(ntim),
      step_factor((ntim - 1)/2.0),
      rho(radius),
      tolerance(toler),
      timestep(2.0/(ntim - 1)),
      future_time(1 + timestep),
      lambdas(nlam),
      times(Eigen::VectorXd::LinSpaced(ntim, -1, 1))
{
  auto xs(linspace(0, M_PI + 2, n_lambda + 1));
  for(int i = 0; i < nlam; ++i) {
    lambdas[i] = rho*semidisk(xs.at(i));
  }

  Eigen::VectorXd predictors(compute_coefficients(predictor_matrix()));
  Eigen::VectorXd correctors(compute_coefficients(corrector_matrix()));

  // Eigen defaults to column major, hence the (n_time x 2) ordering
  predictor_coefs = Eigen::Map<Eigen::ArrayXXd>(predictors.data(), n_time, 2);
  predictor_coefs.col(1) *= step_factor;

  corrector_coefs = Eigen::Map<Eigen::ArrayXXd>(correctors.data(), n_time, 2);
  corrector_coefs.col(1) *= step_factor;
  future_coef = correctors(2*n_time)*step_factor;
}

Eigen::MatrixXcd PredictorCorrector::predictor_matrix() const
{
  Eigen::MatrixXcd result(n_lambda, 2*n_time);

  const Eigen::ArrayXXcd b((lambdas*times.transpose()).array().exp());

  result.block(0, 0, n_lambda, n_time) = b;
  result.block(0, n_time, n_lambda, n_time) = b.colwise()*lambdas.array();

  return result;
}

Eigen::MatrixXcd PredictorCorrector::corrector_matrix() const
{
  Eigen::MatrixXcd result(n_lambda, 2*n_time + 1);

  result.block(0, 0, n_lambda, 2*n_time) = predictor_matrix();
  result.block(0, 2*n_time, n_lambda, 1) = lambdas.array()*rhs_vector().array();

  return result;
}

Eigen::VectorXd PredictorCorrector::compute_coefficients(const Eigen::MatrixXcd &mat) const
{
  Eigen::JacobiSVD<Eigen::MatrixXcd> decomp =
    mat.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).setThreshold(tolerance);

  Eigen::VectorXcd b(rhs_vector()), least_squares(decomp.solve(b));

  return least_squares.real();
}

Eigen::VectorXcd PredictorCorrector::rhs_vector() const
{
  Eigen::ArrayXcd b(lambdas*future_time);
  return b.exp();
}

std::complex<double> semidisk(const double t)
{
  const std::complex<double> iu(0, 1);

  if(0 <= t && t < 1) return iu*t;
  if(1 <= t && t < M_PI + 1) return std::exp(iu*(t + M_PI_2 - 1));
  if(M_PI + 1 <= t && t < M_PI + 2) return iu*(t - (M_PI + 2));

  return 0;
}
