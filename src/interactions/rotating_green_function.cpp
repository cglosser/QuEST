#include "rotating_green_function.h"

std::vector<Eigen::Matrix3cd> GreenFunction::RotatingDyadic::coefficients(
    const Eigen::Vector3d &dr,
    const Interpolation::UniformLagrangeSet &interp) const
{
  std::vector<Eigen::Matrix3cd> coefficients(interp.order() + 1,
                                      Eigen::Matrix3cd::Zero());

  const auto dyads(spatial_dyads(dr));

  for(int i = 0; i <= interp.order(); ++i) {
    coefficients[i] =
        -mu0_over_4pi_hbar_ * std::exp(-iu * omega_ * dr.norm() / c_) *
        (dyads[0].cast<cmplx>() * interp.evaluations[0][i] +
         dyads[1].cast<cmplx>() * (interp.evaluations[1][i] +
                                   iu * omega_ * interp.evaluations[0][i]) +
         dyads[2].cast<cmplx>() *
             (interp.evaluations[2][i] +
              2.0 * iu * omega_ * interp.evaluations[1][i] -
              std::pow(omega_, 2) * interp.evaluations[0][i]));
  }

  return coefficients;
}
