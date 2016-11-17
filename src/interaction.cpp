#include "interaction.h"

Interaction::Interaction(const QuantumDot &d1, const QuantumDot &d2)
{
  Eigen::Vector3d dr(d2.pos - d1.pos);

  const double dimensionless_delay = dr.norm()/(config.c0*config.dt);
  delay = compute_delay(dimensionless_delay);

}

std::pair<int, double> Interaction::compute_delay(const double delay) const
{
  std::pair<int, double> result;

  double idelay;
  result.second = std::modf(delay, &idelay);
  result.first = static_cast<int>(idelay);

  return result;
}


//std::pair<Eigen::Vector3d, Eigen::Vector3d>
  //Interaction::operator()(const QuantumDot &d1, const QuantumDot &d2) const
//{
  //Eigen::Vector3d field1(Eigen::Vector3d::Zero()),
                  //field2(Eigen::Vector3d::Zero());

  //const size_t start = (d1.history.size() - 1) - delay.first;

  //for(int i = 0; i <= config.interpolation_order; ++i) {
    //field1 += farfield_dyadic(d2.polarization(start - i))*interp.weights[2][i] +
              //midfield_dyadic(d2.polarization(start - i))*interp.weights[1][i] +
              //nearfield_dyadic(d2.polarization(start - i))*interp.weights[0][i];

    //field2 += farfield_dyadic(d1.polarization(start - i))*interp.weights[2][i] +
              //midfield_dyadic(d1.polarization(start - i))*interp.weights[1][i] +
              //nearfield_dyadic(d1.polarization(start - i))*interp.weights[0][i];
  //};

  //return std::pair<Eigen::Vector3d, Eigen::Vector3d>(field1, field2);
//}

Eigen::Matrix3d Interaction::rhat_dyadic(const Eigen::Vector3d &dr) const
{
  return dr*dr.transpose()/dr.squaredNorm();
}

//Eigen::Vector3d Interaction::farfield_dyadic(const Eigen::Vector3d &polarization) const
//{
  //return (Eigen::Matrix3d::Identity() - rhat_dyadic())*polarization/dr.norm();
//}

//Eigen::Vector3d Interaction::midfield_dyadic(const Eigen::Vector3d &polarization) const
//{
  //return (Eigen::Matrix3d::Identity() - 3*rhat_dyadic())*polarization*config.c0/dr.squaredNorm();
//}

//Eigen::Vector3d Interaction::nearfield_dyadic(const Eigen::Vector3d &polarization) const
//{
  //double r_cubed = std::pow(dr.norm(), 3);
  //return (Eigen::Matrix3d::Identity() - 3*rhat_dyadic())*polarization*std::pow(config.c0, 2)/r_cubed;
//}
