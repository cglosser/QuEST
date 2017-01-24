#ifndef GREEN_FUNCTION_H
#define GREEN_FUNCTION_H

#include <Eigen/Dense>
#include <array>
#include <vector>

#include "../common.h"
#include "../history.h"
#include "../lagrange_set.h"

namespace GreenFunction {
  class Dyadic;
}

class GreenFunction::Dyadic {
 public:
  Dyadic(const double mu0, const double c, const double hbar)
      : c_(c), mu0_over_4pi_hbar_(mu0 / (4 * M_PI * hbar)){};
  virtual std::vector<Eigen::Matrix3cd> coefficients(
      const Eigen::Vector3d &, const Interpolation::UniformLagrangeSet &) const;

  virtual cmplx polarization_prefactor(
      const History::soltype &matrix_elements) const
  {
    return 2 * matrix_elements[1].real();
  }

 protected:
  double c_, mu0_over_4pi_hbar_;

  std::array<Eigen::Matrix3d, 3> spatial_dyads(const Eigen::Vector3d &dr) const
  {
    std::array<Eigen::Matrix3d, 3> results = {
        {identity_minus_3rsq(dr) * std::pow(c_, 2) / std::pow(dr.norm(), 3),
         identity_minus_3rsq(dr) * c_ / dr.squaredNorm(),
         identity_minus_rsq(dr) / dr.norm()}};

    return results;
  }

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
};

#endif
