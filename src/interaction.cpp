#include "interaction.h"

Interaction::Interaction(const QuantumDot &d1, const QuantumDot &d2)
    : coefs(config.interpolation_order + 1)
{
  const Eigen::Vector3d dr(d2.pos - d1.pos);

  const double dimensionless_delay = dr.norm() / (config.c0 * config.dt);
  delay = compute_delay(dimensionless_delay);

  UniformLagrangeSet interp(config.interpolation_order, delay.second);


  for(size_t i = 0; i < coefs.size(); ++i) {
    coefs.at(i) =
        nearfield_dyadic(dr, d1.dipole, d2.dipole) * interp.weights[0][i] +
        midfield_dyadic(dr, d1.dipole, d2.dipole) * interp.weights[1][i] +
        farfield_dyadic(dr, d1.dipole, d2.dipole) * interp.weights[2][i];
  }
}

std::pair<int, double> Interaction::compute_delay(const double delay) const
{
  std::pair<int, double> result;

  double idelay;
  result.second = std::modf(delay, &idelay);
  result.first = static_cast<int>(idelay);

  return result;
}

Eigen::Matrix3d Interaction::rhat_dyadic(const Eigen::Vector3d &dr) const
{
  return dr * dr.transpose() / dr.squaredNorm();
}

double Interaction::nearfield_dyadic(const Eigen::Vector3d &dr,
                                     const Eigen::Vector3d &src,
                                     const Eigen::Vector3d &obs) const
{
  const double dyad =
    obs.transpose() * (Eigen::Matrix3d::Identity() - 3*rhat_dyadic(dr)) * src;
  return dyad * std::pow(config.c0, 2) / std::pow(dr.norm(), 3);
}

double Interaction::midfield_dyadic(const Eigen::Vector3d &dr,
                                    const Eigen::Vector3d &src,
                                    const Eigen::Vector3d &obs) const
{
  const double dyad =
    obs.transpose() * (Eigen::Matrix3d::Identity() - 3*rhat_dyadic(dr)) * src;
  return dyad * config.c0 / dr.squaredNorm();
}

double Interaction::farfield_dyadic(const Eigen::Vector3d &dr,
                                    const Eigen::Vector3d &src,
                                    const Eigen::Vector3d &obs) const
{
  const double dyad =
    obs.transpose() * (Eigen::Matrix3d::Identity() - rhat_dyadic(dr)) * src;
  return dyad / dr.norm();
}
