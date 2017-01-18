#ifndef GREEN_FUNCTION_H
#define GREEN_FUNCTION_H

#include <Eigen/Dense>
#include <vector>

#include "../lagrange_set.h"

namespace GreenFunction {
  class Dyadic;
}

class GreenFunction::Dyadic {
 public:
  Dyadic(const double mu0, const double c) : mu0_(mu0), c_(c){};
  virtual std::vector<Eigen::Matrix3d> coefficients(
      const Eigen::Vector3d &, const Interpolation::UniformLagrangeSet &) const;

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
    return Eigen::Matrix3d::Identity() - 3 * rhat_dyadic(dr);
  }

 private:
  double mu0_, c_;
};

#endif
