#ifndef GREEN_FUNCTION_H
#define GREEN_FUNCTION_H

#include <Eigen/Dense>
#include <array>
#include <vector>

#include "../common.h"
#include "../lagrange_set.h"

namespace GreenFunction {
  typedef Eigen::Vector2cd soltype;
  class Dyadic;
}

class GreenFunction::Dyadic {
 public:
  Dyadic(const double mu0, const double c, const double hbar)
      : c_(c), mu0_over_4pi_hbar_(mu0 / (4 * M_PI * hbar)){};
  virtual std::vector<Eigen::Matrix3cd> coefficients(
      const Eigen::Vector3d &, const Interpolation::UniformLagrangeSet &) const;

  virtual cmplx polarization_prefactor(const soltype &matrix_elements) const
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

namespace Propagation {
  template <class Derived>
  class Propagator;

  class FixedFramePropagator;
  class RotatingFramePropagator;
}

template <class Derived>
class Propagation::Propagator {
 public:
  Propagator(const double mu0, const double c, const double hbar)
      : c_(c), mu0_over_4pi_hbar_(mu0 / (4 * M_PI * hbar)){};
  auto coefficients(const Eigen::Vector3d &dr,
                    const Interpolation::UniformLagrangeSet &interp) const
  {
    return static_cast<Derived *>(this)->coefficients(dr, interp);
  };

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

class Propagation::FixedFramePropagator
    : public Propagation::Propagator<Propagation::FixedFramePropagator> {
 public:
  FixedFramePropagator(const double mu0, const double c, const double hbar)
      : Propagator(mu0, c, hbar){};
  std::vector<Eigen::Matrix3d> coefficients(
      const Eigen::Vector3d &dr,
      const Interpolation::UniformLagrangeSet &interp) const
  {
    std::vector<Eigen::Matrix3d> coefs(interp.order() + 1,
                                       Eigen::Matrix3d::Zero());

    const auto dyads(spatial_dyads(dr));

    for(int i = 0; i <= interp.order(); ++i) {
      coefs[i] += mu0_over_4pi_hbar_ * (dyads[0] * interp.weights[0][i] +
                                        dyads[1] * interp.weights[1][i] +
                                        dyads[2] * interp.weights[2][i]);
    }

    return coefs;
  }
};

class Propagation::RotatingFramePropagator
    : public Propagation::Propagator<Propagation::RotatingFramePropagator> {
 public:
  RotatingFramePropagator(const double mu0, const double c, const double hbar,
                          const double omega)
      : Propagator(mu0, c, hbar), omega(omega){};
  std::vector<Eigen::Matrix3cd> coefficients(
      const Eigen::Vector3d &dr,
      const Interpolation::UniformLagrangeSet &interp) const
  {
    std::vector<Eigen::Matrix3cd> coefs(interp.order() + 1,
                                        Eigen::Matrix3cd::Zero());

    const auto dyads(spatial_dyads(dr));

    for(int i = 0; i <= interp.order(); ++i) {
      coefs[i] =
          -mu0_over_4pi_hbar_ * std::exp(-iu * omega * dr.norm() / c_) *
          (dyads[0].cast<cmplx>() * interp.weights[0][i] +
           dyads[1].cast<cmplx>() *
               (interp.weights[1][i] + iu * omega * interp.weights[0][i]) +
           dyads[2].cast<cmplx>() *
               (interp.weights[2][i] + 2.0 * iu * omega * interp.weights[1][i] -
                std::pow(omega, 2) * interp.weights[0][i]));
    }

    return coefs;
  }

 private:
  double omega;
};

#endif
