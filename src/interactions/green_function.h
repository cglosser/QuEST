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
  Propagator(const double k2, const double c) : k2_(k2), c_(c){};
  auto coefficients(const Eigen::Vector3d &dr,
                    const Interpolation::UniformLagrangeSet &interp) const
  {
    return static_cast<Derived *>(this)->coefficients(dr, interp);
  };

 protected:
  double k2_, c_;

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

class Propagation::FixedFramePropagator
    : public Propagation::Propagator<Propagation::FixedFramePropagator> {
 public:
  FixedFramePropagator(const double k2, const double c) : Propagator(k2, c){};
  std::vector<Eigen::Matrix3d> coefficients(
      const Eigen::Vector3d &dr,
      const Interpolation::UniformLagrangeSet &interp) const
  {
    std::vector<Eigen::Matrix3d> coefs(interp.order() + 1,
                                       Eigen::Matrix3d::Zero());

    const auto dyads(spatial_dyads(dr));

    for(int i = 0; i <= interp.order(); ++i) {
      coefs[i] = -k2_ * (dyads[0] * interp.evaluations[0][i] +
                         dyads[1] * interp.evaluations[1][i] +
                         dyads[2] * interp.evaluations[2][i]);
    }

    return coefs;
  }
};

class Propagation::RotatingFramePropagator
    : public Propagation::Propagator<Propagation::RotatingFramePropagator> {
 public:
  RotatingFramePropagator(const double k2, const double c, const double omega)
      : Propagator(k2, c), omega(omega){};
  std::vector<Eigen::Matrix3cd> coefficients(
      const Eigen::Vector3d &dr,
      const Interpolation::UniformLagrangeSet &interp) const
  {
    std::vector<Eigen::Matrix3cd> coefs(interp.order() + 1,
                                        Eigen::Matrix3cd::Zero());

    const auto dyads(spatial_dyads(dr));

    for(int i = 0; i <= interp.order(); ++i) {
      coefs[i] =
          -k2_ * std::exp(-iu * omega * dr.norm() / c_) *
          (dyads[0].cast<cmplx>() * interp.evaluations[0][i] +
           dyads[1].cast<cmplx>() * (interp.evaluations[1][i] +
                                     iu * omega * interp.evaluations[0][i]) +
           dyads[2].cast<cmplx>() *
               (interp.evaluations[2][i] +
                2.0 * iu * omega * interp.evaluations[1][i] -
                std::pow(omega, 2) * interp.evaluations[0][i]));
    }

    return coefs;
  }

 private:
  double omega;
};

#endif
