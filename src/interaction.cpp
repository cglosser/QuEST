#include "interaction.h"

Interaction::Interaction(const Eigen::Vector3d &sep) : dr(sep)
{
  double dimensionless_delay = dr.norm()/(config.c0*config.dt);

  double idelay;
  delay.second = std::modf(dimensionless_delay, &idelay);
  delay.first = static_cast<int>(idelay);

  interp = UniformLagrangeSet(delay.second);
}

std::pair<Eigen::Vector3d, Eigen::Vector3d>
  Interaction::operator()(const QuantumDot &d1, const QuantumDot &d2) const
{
  Eigen::Vector3d field1(Eigen::Vector3d::Zero()),
                  field2(Eigen::Vector3d::Zero());

  return std::pair<Eigen::Vector3d, Eigen::Vector3d>(field1, field2);
}

Eigen::Matrix3d Interaction::rhat_dyadic() const
{
  return dr*dr.transpose()/dr.squaredNorm();
}
