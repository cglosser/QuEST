#ifndef RHS_H
#define RHS_H

#include <memory>
#include "../history.h"

namespace Integrator {
  template <class soltype>
  class RHS;
}

template <class soltype>
class Integrator::RHS {
 public:
  RHS(const double dt, const std::shared_ptr<History<soltype>> &history)
      : dt(dt), history(history){};
  virtual void evaluate(const int) const = 0;

 protected:
  double dt;
  std::shared_ptr<History<soltype>> history;
};

#endif
