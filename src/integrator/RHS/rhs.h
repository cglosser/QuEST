#ifndef RHS_H
#define RHS_H

#include <functional>
#include <memory>
#include <vector>
#include "../history.h"

namespace Integrator {
  template <class soltype>
  class RHS;
}

template <class soltype>
class Integrator::RHS {
 public:
  RHS(const double dt, const std::shared_ptr<History<soltype>> &history,
      const std::vector<std::function<soltype(soltype, double)>> &rhs_functions)
      : dt(dt), history(history), rhs_functions(rhs_functions){};
  void evaluate(const int) const;

 protected:
  double dt;
  std::shared_ptr<History<soltype>> history;
  std::vector<std::function<double(soltype, double)>> rhs_functions;
};

template <class soltype>
void Integrator::RHS<soltype>::evaluate(const int n) const
{
  const double time = n * dt;
  for(int i = 0; i < static_cast<int>(history->array.shape()[0]); ++i) {
    history->array[i][n][1] = rhs_functions[i](history->array[i][n][0], time);
  }
}

#endif
