#include "integrator.h"

std::complex<double> semidisk(const double t)
{
  const std::complex<double> iu(0, 1);

  if(0 <= t && t < 1) return iu * t;
  if(1 <= t && t < M_PI + 1) return std::exp(iu * (t + M_PI_2 - 1));
  if(M_PI + 1 <= t && t < M_PI + 2) return iu * (t - (M_PI + 2));

  return 0;
}

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

PredictorCorrector::Weights::Weights(const size_t n_lambda, const size_t n_time,
                                     const double radius,
                                     const double tolerance)
    : n_time_(n_time)
{
  const double step_factor = (n_time - 1) / 2.0;
  WeightsBuilder builder(n_lambda, n_time, radius, tolerance);

  Eigen::VectorXd predictors(builder.predictors());
  Eigen::VectorXd correctors(builder.correctors());

  // Eigen defaults to column major, hence the (n_time x 2) ordering
  ps = Eigen::Map<Eigen::ArrayXXd>(predictors.data(), n_time, 2).transpose();
  cs = Eigen::Map<Eigen::ArrayXXd>(correctors.data(), n_time, 2).transpose();

  ps.row(1) *= step_factor;
  cs.row(1) *= step_factor;
  future_coef = correctors(2 * n_time) * step_factor;
}

PredictorCorrector::Integrator::Integrator(
    std::vector<QuantumDot> qds, const size_t n, const double timestep,
    const size_t n_lambda, const size_t n_time, const double radius,
    const double tolerance)
    : num_steps(n + 1),
      dt(timestep),
      weights(n_lambda, n_time, radius, tolerance),
      history(
          boost::extents[qds.size()]
                        [HistoryArray::extent_range(-n_time, num_steps)][2]),
      rabi_freqs(qds.size(), 0)
{
  dots.swap(qds);

  for(int dot_idx = 0; dot_idx < static_cast<int>(dots.size()); ++dot_idx) {
    for(int i = -weights.width(); i <= 0; ++i) {
      history[dot_idx][i][0] = matrix_elements(1, 0);
      history[dot_idx][i][1] = matrix_elements(0, 0);
    }
  }

  now = 1;
}

void PredictorCorrector::Integrator::step()
{
  assert(now < static_cast<int>(num_steps));

  predictor();

  set_rabi_freqs(now * dt);
  evaluator();

  for(int m = 1; m < 10; ++m) {
    corrector();

    set_rabi_freqs(now * dt);
    evaluator();
  }

  ++now;
}

void PredictorCorrector::Integrator::predictor()
{
  const int start = now - weights.width();

  for(int sol_idx = 0; sol_idx < static_cast<int>(dots.size()); ++sol_idx) {
    history[sol_idx][now][0].setZero();
    for(int h = 0; h < static_cast<int>(weights.width()); ++h) {
      history[sol_idx][now][0] +=
          history[sol_idx][start + h][0] * weights.ps(0, h) +
          history[sol_idx][start + h][1] * weights.ps(1, h) * dt;
    }
  }
}

void PredictorCorrector::Integrator::corrector()
{
  const int start = now - weights.width();

  for(int sol_idx = 0; sol_idx < static_cast<int>(dots.size()); ++sol_idx) {
    history[sol_idx][now][0] =
          weights.future_coef * history[sol_idx][now][1] * dt;
    for(int h = 0; h < static_cast<int>(weights.width()); ++h) {
      history[sol_idx][now][0] +=
          history[sol_idx][start + h][0] * weights.cs(0, h) +
          history[sol_idx][start + h][1] * weights.cs(1, h) * dt;
    }
  }
}

void PredictorCorrector::Integrator::evaluator()
{
  for(size_t sol_idx = 0; sol_idx < dots.size(); ++sol_idx) {
    history[sol_idx][now][1] = dots[sol_idx].liouville_rhs(
        history[sol_idx][now][0], rabi_freqs[sol_idx]);
  }
}

void PredictorCorrector::Integrator::set_rabi_freqs(const double time)
{
  const double chi_0 = M_PI / std::sqrt(2 * M_PI * std::pow(0.1, 2));
  const double x = chi_0 * gaussian((time - 0.5) / 0.1) * cos(2278.9 * time);
  std::fill(rabi_freqs.begin(), rabi_freqs.end(), x);
}
