#include "fixed_dyadic.h"

GreenFunction::FixedFrameDyadic::FixedFrameDyadic(const double mu0, const double c) : mu0_(mu0), c_(c)
{
}

std::vector<Eigen::Matrix3cd>
GreenFunction::FixedFrameDyadic::calculate_coefficients(
    const Eigen::Vector3d &dr, const Interpolation::UniformLagrangeSet &interp) const
{
  std::vector<Eigen::Matrix3cd> coefficients(interp.order() + 1);

  const double alpha = -mu0_ / (4 * M_PI);
  const Eigen::Matrix3d nearfield(identity_minus_3rsq(dr) * std::pow(c_, 2) /
                                  std::pow(dr.norm(), 3));
  const Eigen::Matrix3d midfield(identity_minus_3rsq(dr) * c_ /
                                 dr.squaredNorm());
  const Eigen::Matrix3d farfield(identity_minus_rsq(dr) / dr.norm());

  for(int i = 0; i <= interp.order(); ++i) {
    coefficients[i] =
        alpha *
        (nearfield * interp.weights[0][i] + midfield * interp.weights[1][i] +
         farfield * interp.weights[2][i])
            .cast<std::complex<double>>();
  }

  return coefficients;
}
