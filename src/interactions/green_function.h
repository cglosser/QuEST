#ifndef GREEN_FUNCTION_H
#define GREEN_FUNCTION_H

#include <Eigen/Dense>
#include <vector>

#include "../lagrange_set.h"

namespace GreenFunction {
  class TimeDomainGreenFunction;
}

class GreenFunction::TimeDomainGreenFunction{
 public:
  virtual std::vector<Eigen::Matrix3cd> calculate_coefficients(
      const Eigen::Vector3d &, const Interpolation::UniformLagrangeSet &) const = 0;

  static Eigen::Matrix3d rhat_dyadic(const Eigen::Vector3d &dr)
  {
    return dr * dr.transpose() / dr.squaredNorm();
  }

  static Eigen::Matrix3d identity_minus_rsq(const Eigen::Vector3d &dr)
  {
    return Eigen::Matrix3d::Identity() - rhat_dyadic(dr);
  }

  static Eigen::Matrix3d identity_minus_3rsq(const Eigen::Vector3d &dr)
  {
    return Eigen::Matrix3d::Identity() - 3*rhat_dyadic(dr);
  }
};

#endif
