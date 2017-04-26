#include "integrator.h"

constexpr int NUM_CORRECTOR_STEPS = 10;

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
  WeightsBuilder(const int, const int, const double);

  Eigen::VectorXd predictors() { return compute_coeff(predictor_matrix()); }
  Eigen::VectorXd correctors() { return compute_coeff(corrector_matrix()); }
 private:
  Eigen::VectorXcd lambdas;
  Eigen::VectorXd times;
  double timestep, future_time;

  Eigen::MatrixXcd predictor_matrix() const;
  Eigen::MatrixXcd corrector_matrix() const;
  Eigen::VectorXcd rhs_vector() const;

  Eigen::VectorXd compute_coeff(const Eigen::MatrixXcd &) const;
};

WeightsBuilder::WeightsBuilder(const int n_lambda, const int n_time,
                               const double radius)
    : lambdas(n_lambda),
      times(Eigen::VectorXd::LinSpaced(n_time, -1, 1)),
      timestep(2.0 / (n_time - 1)),
      future_time(1 + timestep)
{
  Eigen::VectorXd xs(Eigen::VectorXd::LinSpaced(n_lambda + 1, 0, M_PI + 2));
  for(int i = 0; i < n_lambda; ++i) lambdas[i] = radius * semidisk(xs[i]);
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
      mat.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV);

  Eigen::VectorXcd b(rhs_vector()), least_squares(decomp.solve(b));

  return least_squares.real();
}

PredictorCorrector::Weights::Weights(const int n_lambda, const int n_time,
                                     const double radius)
    : n_time(n_time)
{
  const double step_factor = (n_time - 1) / 2.0;
  WeightsBuilder builder(n_lambda, n_time, radius);

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
    const double dt, const int n_lambda, const int n_time, const double radius,
    const std::shared_ptr<History::HistoryArray> &history,
    const std::vector<rhs_func> &rhs_funcs,
    std::vector<std::shared_ptr<Interaction>> interactions)
    : num_solutions(rhs_funcs.size()),
      time_idx_ubound(history->index_bases()[1] + history->shape()[1]),
      dt(dt),
      weights(n_lambda, n_time, radius),
      history(history),
      rhs_funcs(rhs_funcs),
      interactions(std::move(interactions))
{
  assert(rhs_funcs.size() == history->shape()[0]);
}

void PredictorCorrector::Integrator::solve() const
{
  for(int step = 0; step < time_idx_ubound; ++step) {
    solve_step(step);
    log_percentage_complete(step);
    throw_if_unbounded_solution(step);
  }
}

void PredictorCorrector::Integrator::solve_step(const int step) const
{
  assert(0 <= step && step < time_idx_ubound);

  predictor(step);
  evaluator(step);

  for(int m = 0; m < NUM_CORRECTOR_STEPS; ++m) {
    corrector(step);
    evaluator(step);
  }
}

void PredictorCorrector::Integrator::predictor(const int step) const
{
  const int start = step - weights.width();

  for(int sol_idx = 0; sol_idx < num_solutions; ++sol_idx) {
    (*history)[sol_idx][step][0].setZero();
    for(int h = 0; h < static_cast<int>(weights.width()); ++h) {
      (*history)[sol_idx][step][0] +=
          (*history)[sol_idx][start + h][0] * weights.ps(0, h) +
          (*history)[sol_idx][start + h][1] * weights.ps(1, h) * dt;
    }
  }
}

void PredictorCorrector::Integrator::corrector(const int step) const
{
  const int start = step - weights.width();

  for(int sol_idx = 0; sol_idx < num_solutions; ++sol_idx) {
    (*history)[sol_idx][step][0] =
        weights.future_coef * (*history)[sol_idx][step][1] * dt;
    for(int h = 0; h < static_cast<int>(weights.width()); ++h) {
      (*history)[sol_idx][step][0] +=
          (*history)[sol_idx][start + h][0] * weights.cs(0, h) +
          (*history)[sol_idx][start + h][1] * weights.cs(1, h) * dt;
    }
  }
}

void PredictorCorrector::Integrator::evaluator(const int step) const
{
  auto projected_efields =
      std::accumulate(interactions.begin(), interactions.end(),
                      Interaction::ResultArray::Zero(num_solutions, 1).eval(),
                      [step](const Interaction::ResultArray &r,
                             const std::shared_ptr<Interaction> &interaction) {
                        return r + interaction->evaluate(step);
                      });

  for(int solution = 0; solution < num_solutions; ++solution) {
    (*history)[solution][step][1] = rhs_funcs[solution](
        (*history)[solution][step][0], projected_efields[solution]);
  }
}

bool PredictorCorrector::Integrator::all_finite(const int step) const
{
  for(int sol_idx = 0; sol_idx < num_solutions; ++sol_idx) {
    History::soltype &sol = (*history)[sol_idx][step][0];
    History::soltype &dsol = (*history)[sol_idx][step][1];

    if(!History::isfinite(sol) || !History::isfinite(dsol)) {
      return false;
    }
  }
  return true;
}

void PredictorCorrector::Integrator::throw_if_unbounded_solution(
    const int step) const
{
  if(!all_finite(step)) {
    const std::string msg = "unbounded history values at or before step ";
    throw std::domain_error(msg + std::to_string(step));
  }
}

void PredictorCorrector::Integrator::log_percentage_complete(
    const int step) const
{
  if(step % (time_idx_ubound / 10) == 0) {
    std::cout << "\t" << static_cast<double>(step) / time_idx_ubound
              << std::endl;
  }
}
