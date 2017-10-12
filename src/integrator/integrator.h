#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <string>

#include "../math_utils.h"
#include "RHS/rhs.h"
#include "history.h"
#include "weights.h"

namespace Integrator {
  template <class soltype>
  class PredictorCorrector;
}

template <class soltype>
class Integrator::PredictorCorrector {
 public:
  PredictorCorrector(const double, const int, const int, const double,
                     const std::shared_ptr<Integrator::History<soltype>> &,
                     std::unique_ptr<Integrator::RHS<soltype>>);
  void solve() const;

 private:
  int num_solutions, time_idx_ubound;
  double dt;
  Weights weights;
  std::shared_ptr<Integrator::History<soltype>> history;
  std::unique_ptr<Integrator::RHS<soltype>> rhs;

  void solve_step(const int) const;
  void predictor(const int) const;
  void corrector(const int) const;
};

constexpr int NUM_CORRECTOR_STEPS = 10;

template <class soltype>
Integrator::PredictorCorrector<soltype>::PredictorCorrector(
    const double dt, const int n_lambda, const int n_time, const double radius,
    const std::shared_ptr<Integrator::History<soltype>> &history,
    std::unique_ptr<Integrator::RHS<soltype>> rhs)
    : num_solutions(history->array.shape()[0]),
      time_idx_ubound(history->array.index_bases()[1] +
                      history->array.shape()[1]),
      dt(dt),
      weights(n_lambda, n_time, radius),
      history(history),
      rhs(std::move(rhs))
{
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::solve() const
{
  for(int step = 0; step < time_idx_ubound; ++step) {
    solve_step(step);
  }
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::solve_step(const int step) const
{
  assert(0 <= step && step < time_idx_ubound);

  predictor(step);
  rhs->evaluate(step);

  for(int m = 0; m < NUM_CORRECTOR_STEPS; ++m) {
    corrector(step);
    rhs->evaluate(step);
  }
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::predictor(const int step) const
{
  const int start = step - weights.width();

  for(int sol_idx = 0; sol_idx < num_solutions; ++sol_idx) {
    for(int h = 0; h < static_cast<int>(weights.width()); ++h) {
      history->array[sol_idx][step][0] +=
          history->array[sol_idx][start + h][0] * weights.ps(0, h) +
          history->array[sol_idx][start + h][1] * weights.ps(1, h) * dt;
    }
  }
}

template <class soltype>
void Integrator::PredictorCorrector<soltype>::corrector(const int step) const
{
  const int start = step - weights.width();

  for(int sol_idx = 0; sol_idx < num_solutions; ++sol_idx) {
    history->array[sol_idx][step][0] =
        weights.future_coef * history->array[sol_idx][step][1] * dt;
    for(int h = 0; h < static_cast<int>(weights.width()); ++h) {
      history->array[sol_idx][step][0] +=
          history->array[sol_idx][start + h][0] * weights.cs(0, h) +
          history->array[sol_idx][start + h][1] * weights.cs(1, h) * dt;
    }
  }
}

#endif
