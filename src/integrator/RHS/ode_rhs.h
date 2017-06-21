#ifndef ODE_RHS_H
#define ODE_RHS_H

#include <cmath>
#include "rhs.h"

namespace Integrator {
  class ODE_RHS;
}

class Integrator::ODE_RHS : public RHS<double> {
 public:
  ODE_RHS(const double dt,
          const std::shared_ptr<Integrator::History<double>> &history)
      : RHS(dt, history){};
  void evaluate(const int) const;

 private:
};

void Integrator::ODE_RHS::evaluate(const int n) const
{
  const double time = n * dt;
  for(int i = 0; i < static_cast<int>(history->array.shape()[0]); ++i) {
    history->array[i][n][1] =
        1 / (1 + std::exp(-(time - 10))) + history->array[i][n][0];
  }
}

#endif
