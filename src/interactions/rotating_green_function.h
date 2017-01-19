#ifndef ROTATING_GREEN_FUNCTION_H
#define ROTATING_GREEN_FUNCTION_H

#include "green_function.h"

namespace GreenFunction {
  class RotatingDyadic;
}

class GreenFunction::RotatingDyadic : public GreenFunction::Dyadic {
 public:
  RotatingDyadic(const double mu0, const double c, const double omega)
      : Dyadic(mu0, c), omega_(omega){};
  std::vector<Eigen::Matrix3cd> coefficients(
      const Eigen::Vector3d &, const Interpolation::UniformLagrangeSet &) const;

 private:
  double omega_;
};

#endif
