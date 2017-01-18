#include "green_function.h"

std::vector<Eigen::Matrix3cd> GreenFunction::Dyadic::coefficients(
    const Eigen::Vector3d &dr,
    const Interpolation::UniformLagrangeSet &interp) const
{
  std::vector<Eigen::Matrix3cd> coefs(interp.order() + 1);

  const double alpha = -mu0_ / (4 * M_PI);
  const Eigen::Matrix3d nearfield(identity_minus_3rsq(dr) * std::pow(c_, 2) /
                                  std::pow(dr.norm(), 3));
  const Eigen::Matrix3d midfield(identity_minus_3rsq(dr) * c_ /
                                 dr.squaredNorm());
  const Eigen::Matrix3d farfield(identity_minus_rsq(dr) / dr.norm());

  for(int i = 0; i <= interp.order(); ++i) {
    Eigen::Matrix3d dyad = alpha * (nearfield * interp.weights[0][i] +
                                    midfield * interp.weights[1][i] +
                                    farfield * interp.weights[2][i]);
    coefs[i] = dyad.cast<std::complex<double>>();
  }

  return coefs;
}
