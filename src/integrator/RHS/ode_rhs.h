#ifndef ODE_RHS_H
#define ODE_RHS_H

#include <cmath>
#include "rhs.h"

namespace Integrator {
  class ODE_RHS;
}

class Integrator::ODE_RHS : public RHS<double> {
 public:
  ODE_RHS(double dt) : RHS<double>(dt){};
  double operator()(const double &, const int) const;

 private:
};

double Integrator::ODE_RHS::operator()(const double &f, const int n) const
{
  const double time = n * dt;
  return 1 / (1 + std::exp(-(time - 10))) + f;
}

#endif
