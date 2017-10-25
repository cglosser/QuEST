#include "euler.h"

EulerIntegrator::EulerIntegrator(
    const int step,
    const double dt,
    const std::shared_ptr<Integrator::History<Eigen::Vector3d>> &history,
    std::unique_ptr<Integrator::RHS<Eigen::Vector3d>> &rhs_functions)
    : step(step), dt(dt), history(history), rhs_functions(rhs_functions)
{
}

void EulerIntegrator::solve() const
{
  rhs_functions->evaluate(step);
  for(unsigned int num = 0; num < history->array.size(); ++num) {
    history->array[step + 1][num][0] =
        history->array[step][num][1] * dt + history->array[step][num][0];
  }
}
