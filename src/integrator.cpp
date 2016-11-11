#include "integrator.h"

using namespace std;

std::complex<double> semidisk(const double);

PredictorCorrector::PredictorCorrector(const int nlam, const int ntim,
  const double radius, const double toler)
    : n_lambda(nlam),
      n_time(ntim),
      rho(radius),
      tolerance(toler),
      timestep(2.0/(ntim - 1)),
      future_time(1 + timestep),
      step_factor((ntim - 1)/2.0),
      lambdas(nlam),
      times(Eigen::VectorXd::LinSpaced(ntim, -1, 1))
{
  auto xs(linspace(0, M_PI + 2, n_lambda + 1));
  for(int i = 0; i < nlam; ++i) {
    lambdas[i] = rho*semidisk(xs.at(i));
  }
}

Eigen::MatrixXcd PredictorCorrector::predictor_matrix() const
{
  Eigen::MatrixXcd result(n_lambda, 2*n_time);

  const Eigen::ArrayXXcd b((lambdas*times.transpose()).array().exp());

  result.block(0, 0, n_lambda, n_time) = b;
  result.block(0, n_time, n_lambda, n_time) = b.colwise()*lambdas.array();

  return result;
}

std::complex<double> semidisk(const double t)
{
  const std::complex<double> iu(0, 1);

  if(0 <= t && t < 1) return iu*t;
  if(1 <= t && t < M_PI + 1) return std::exp(iu*(t + M_PI_2 - 1));
  if(M_PI + 1 <= t && t < M_PI + 2) return iu*(t - (M_PI + 2));

  return 0;
}
