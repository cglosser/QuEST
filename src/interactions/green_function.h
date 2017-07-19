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
  Propagator(const double e0, const double c, const double hbar)
      : c(c), e0(e0), hbar(hbar){};
  auto coefficients(const Eigen::Vector3d &dr,
                    const Interpolation::UniformLagrangeSet &interp) const
  {
    return static_cast<Derived *>(this)->coefficients(dr, interp);
  };

 protected:
  double c, e0, hbar;
};

class Propagation::FixedFramePropagator
    : public Propagation::Propagator<Propagation::FixedFramePropagator> {
 public:
  FixedFramePropagator(const double e0, const double c, const double hbar)
      : Propagator(e0, c, hbar){};
  std::vector<Eigen::Matrix3d> coefficients(
      const Eigen::Vector3d &dr,
      const Interpolation::UniformLagrangeSet &interp) const
  {
    std::vector<Eigen::Matrix3d> coefs(interp.order() + 1,
                                       Eigen::Matrix3d::Zero());

    Eigen::Matrix3d r_matrix;
    r_matrix << 0, dr[2], -dr[1], -dr[2], 0, dr[0], dr[1], -dr[0], 0;

    for(int i = 0; i <= interp.order(); ++i) {
      coefs[i] += e0 / (4 * M_PI) * -r_matrix * 
                  (interp.evaluations[1][i] / std::pow(dr.norm(), 3) -
                   interp.evaluations[2][i] * c / std::pow(dr.norm(), 2));
    }

    return coefs;
  }
};

#endif
