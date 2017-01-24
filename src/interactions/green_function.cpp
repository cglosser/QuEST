#include "green_function.h"

std::vector<Eigen::Matrix3cd> GreenFunction::Dyadic::coefficients(
    const Eigen::Vector3d &dr,
    const Interpolation::UniformLagrangeSet &interp) const
{
  std::vector<Eigen::Matrix3cd> coefs(interp.order() + 1,
                                      Eigen::Matrix3cd::Zero());

  const auto dyads(spatial_dyads(dr));

  for(int i = 0; i <= interp.order(); ++i) {
    for(int term = 0; term < 3; ++term) {
      coefs[i] += mu0_over_4pi_hbar_ * dyads[term].cast<cmplx>() *
                  interp.weights[term][i];
    }
  }

  return coefs;
}
