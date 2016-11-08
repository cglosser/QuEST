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

  const size_t start = d1.history.size() - 1 - delay.first;

  for(int i = 0; i <= config.interpolation_order; ++i) {
    field1 +=    farfield_dyadic(d2.polarization(start - i))*interp.weights[2][config.interpolation_order - i]
              +  midfield_dyadic(d2.polarization(start - i))*interp.weights[1][config.interpolation_order - i]
              + nearfield_dyadic(d2.polarization(start - i))*interp.weights[0][config.interpolation_order - i];

    field2 +=    farfield_dyadic(d1.polarization(start - i))*interp.weights[2][config.interpolation_order - i]
              +  midfield_dyadic(d1.polarization(start - i))*interp.weights[1][config.interpolation_order - i]
              + nearfield_dyadic(d1.polarization(start - i))*interp.weights[0][config.interpolation_order - i];
  };

  return std::pair<Eigen::Vector3d, Eigen::Vector3d>(field1, field2);
}

Eigen::Matrix3d Interaction::rhat_dyadic() const
{
  return dr*dr.transpose()/dr.squaredNorm();
}

Eigen::Vector3d Interaction::farfield_dyadic(const Eigen::Vector3d &polarization) const
{
  return (Eigen::Matrix3d::Identity() - rhat_dyadic())*polarization/dr.norm();
}

Eigen::Vector3d Interaction::midfield_dyadic(const Eigen::Vector3d &polarization) const
{
  return (Eigen::Matrix3d::Identity() - 3*rhat_dyadic())*polarization*config.c0/dr.squaredNorm();
}

Eigen::Vector3d Interaction::nearfield_dyadic(const Eigen::Vector3d &polarization) const
{
  double r_cubed = std::pow(dr.norm(), 3);
  return (Eigen::Matrix3d::Identity() - 3*rhat_dyadic())*polarization*std::pow(config.c0, 2)/r_cubed;
}
