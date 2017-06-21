#ifndef RHS_H
#define RHS_H

#include <memory>

namespace Integrator {
  template <class soltype>
  class RHS;
}

template <class soltype>
class Integrator::RHS {
 public:
  RHS(const double dt, const std::shared_ptr<History::HistoryArray> &history)
      : dt(dt), history(history){};
  virtual void evaluate(const int) = 0;

 protected:
  double dt;
  std::shared_ptr<History::HistoryArray> history;
};

#endif
