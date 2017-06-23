#ifndef ODE_RHS_H
#define ODE_RHS_H

#include <cmath>
#include <functional>
#include "rhs.h"

namespace Integrator {
  class ODE_RHS;
}

class Integrator::ODE_RHS : public RHS<double> {
 public:
  ODE_RHS(const double, const std::shared_ptr<Integrator::History<double>> &,
          const std::function<double(double, double)> &);
  void evaluate(const int) const;

 private:
  std::function<double(double, double)> rhs_func;
};

Integrator::ODE_RHS::ODE_RHS(
    const double dt,
    const std::shared_ptr<Integrator::History<double>> &history,
    const std::function<double(double, double)> &rhs_func)
    : RHS(dt, history), rhs_func(rhs_func)
{
}

void Integrator::ODE_RHS::evaluate(const int n) const
{
  const double time = n * dt;
  for(int i = 0; i < static_cast<int>(history->array.shape()[0]); ++i) {
    history->array[i][n][1] = rhs_func(history->array[i][n][0], time);
  }
}

#endif
