#ifndef WEIGHTS_H
#define WEIGHTS_H

#include <Eigen/Dense>

namespace Integrator{
  class Weights;
}

class Integrator::Weights {
 public:
  Weights(const int, const int, const double);

  Eigen::ArrayXXd ps, cs;
  double future_coef;

  int width() const { return n_time; }
 private:
  int n_time;
};

#endif
