#ifndef RHS_H
#define RHS_H

namespace Integrator {
  template <class soltype>
  class RHS;
}

template <class soltype>
class Integrator::RHS {
 public:
  RHS(const double dt) : dt(dt) {};
  virtual soltype operator()(const soltype &, const int) const = 0;

 protected:
  double dt;
};

#endif
