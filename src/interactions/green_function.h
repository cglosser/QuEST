#ifndef GREEN_FUNCTION_H
#define GREEN_FUNCTION_H

#include <Eigen/Dense>
#include <array>
#include <vector>

#include "../common.h"
#include "../lagrange_set.h"

namespace Propagation {
  template <class Derived>
  class Propagator;

  class FixedFramePropagator;
  class RotatingFramePropagator;
}

template <class Derived>
class Propagation::Propagator {
 public:
  Propagator(const double c) : c(c){};
  auto coefficients(const Eigen::Vector3d &dr,
                    const Interpolation::UniformLagrangeSet &interp) const
  {
    return static_cast<Derived *>(this)->coefficients(dr, interp);
  };

 protected:
  double c;
};

class Propagation::FixedFramePropagator
    : public Propagation::Propagator<Propagation::FixedFramePropagator> {
 public:
  FixedFramePropagator(const double c) : Propagator(c){};
  std::vector<Eigen::Matrix3d> coefficients(
      const Eigen::Vector3d &dr,
      const Interpolation::UniformLagrangeSet &interp) const
  {
    std::vector<Eigen::Matrix3d> coefs(interp.order() + 1,
                                       Eigen::Matrix3d::Zero());

    std::array<Eigen::Matrix3d, 3> dyads(spatial_dyads(dr));

    for(int i = 0; i <= interp.order(); ++i) {
      for(int term = 0; term < 3; ++term) {
        coefs[i] += -1 * dyads[term] * interp.evaluations[term][i];
      }
    }
    return coefs;
  }

 protected:
  std::array<Eigen::Matrix3d, 3> spatial_dyads(const Eigen::Vector3d &dr) const
  {
    std::array<Eigen::Matrix3d, 3> results = {
        {identity_minus_3rsq(dr) / std::pow(dr.norm(), 3),
         identity_minus_3rsq(dr) / (c * dr.squaredNorm()),
         identity_minus_rsq(dr) / (std::pow(c, 2) * dr.norm())}};

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
