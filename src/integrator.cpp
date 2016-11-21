#include "integrator.h"

std::complex<double> semidisk(const double);

class WeightsBuilder {
 public:
  WeightsBuilder(const size_t, const size_t, const double, const double);

  Eigen::VectorXd predictors() { return compute_coeff(predictor_matrix()); }
  Eigen::VectorXd correctors() { return compute_coeff(corrector_matrix()); }
 private:
  Eigen::VectorXcd lambdas;
  Eigen::VectorXd times;
  double timestep, future_time, tolerance;

  Eigen::MatrixXcd predictor_matrix() const;
  Eigen::MatrixXcd corrector_matrix() const;
  Eigen::VectorXcd rhs_vector() const;

  Eigen::VectorXd compute_coeff(const Eigen::MatrixXcd &) const;
};

WeightsBuilder::WeightsBuilder(const size_t n_lambda, const size_t n_time,
                     const double radius, const double toler)
    : lambdas(n_lambda),
      times(Eigen::VectorXd::LinSpaced(n_time, -1, 1)),
      timestep(2.0 / (n_time - 1)),
      future_time(1 + timestep),
      tolerance(toler)
{
  Eigen::VectorXd xs(Eigen::VectorXd::LinSpaced(n_lambda + 1, 0, M_PI + 2));
  for(size_t i = 0; i < n_lambda; ++i) lambdas[i] = radius * semidisk(xs[i]);
}

Eigen::MatrixXcd WeightsBuilder::predictor_matrix() const
{
  Eigen::MatrixXcd result(lambdas.size(), 2 * times.size());

  const Eigen::ArrayXXcd b((lambdas * times.transpose()).array().exp());

  result.block(0, 0, lambdas.size(), times.size()) = b;
  result.block(0, times.size(), lambdas.size(), times.size()) =
      b.colwise() * lambdas.array();

  return result;
}

Eigen::MatrixXcd WeightsBuilder::corrector_matrix() const
{
  Eigen::MatrixXcd result(lambdas.size(), 2 * times.size() + 1);

  result.block(0, 0, lambdas.size(), 2 * times.size()) = predictor_matrix();
  result.block(0, 2 * times.size(), lambdas.size(), 1) =
      lambdas.array() * rhs_vector().array();

  return result;
}

Eigen::VectorXcd WeightsBuilder::rhs_vector() const
{
  Eigen::ArrayXcd b(lambdas * future_time);
  return b.exp();
}

Eigen::VectorXd WeightsBuilder::compute_coeff(const Eigen::MatrixXcd &mat) const
{
  Eigen::JacobiSVD<Eigen::MatrixXcd> decomp =
      mat.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
          .setThreshold(tolerance);

  Eigen::VectorXcd b(rhs_vector()), least_squares(decomp.solve(b));

  return least_squares.real();
}

Weights::Weights(const size_t n_lambda,
                 const size_t n_time, const double radius,
                 const double tolerance)
    : n_time_(n_time)
{
  const double step_factor = (n_time - 1) / 2.0;
  WeightsBuilder builder(n_lambda, n_time, radius, tolerance);

  Eigen::VectorXd predictors(builder.predictors());
  Eigen::VectorXd correctors(builder.correctors());

  // Eigen defaults to column major, hence the (n_time x 2) ordering
  predictor_coefs = Eigen::Map<Eigen::ArrayXXd>(predictors.data(), n_time, 2);
  corrector_coefs = Eigen::Map<Eigen::ArrayXXd>(correctors.data(), n_time, 2);

  predictor_coefs.col(1) *= step_factor;
  corrector_coefs.col(1) *= step_factor;
  future_coef = correctors(2 * n_time) * step_factor;
}

std::complex<double> semidisk(const double t)
{
  const std::complex<double> iu(0, 1);

  if(0 <= t && t < 1) return iu * t;
  if(1 <= t && t < M_PI + 1) return std::exp(iu * (t + M_PI_2 - 1));
  if(M_PI + 1 <= t && t < M_PI + 2) return iu * (t - (M_PI + 2));

  return 0;
}
