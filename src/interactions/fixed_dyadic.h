#ifndef FIXED_DYADIC_H
#define FIXED_DYADIC_H

#include "green_function.h"

namespace GreenFunction {
  class FixedFrameDyadic;
}

class GreenFunction::FixedFrameDyadic
    : public GreenFunction::TimeDomainGreenFunction {
 public:
  FixedFrameDyadic(const double, const double);
  std::vector<Eigen::Matrix3cd> calculate_coefficients(
      const Eigen::Vector3d &, const Interpolation::UniformLagrangeSet &) const;

 private:
  double mu0_, c_;
};

#endif
