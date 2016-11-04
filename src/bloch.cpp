#include "bloch.h"

QuantumDot::QuantumDot(const double omega, const double kappa, 
    const std::pair<double, double> &ts, const Eigen::Vector3d &loc) :
  history(500),
  frequency(omega),
  dipole(kappa),
  damping(ts),
  location(loc)
{}

std::ostream &operator<<(std::ostream &os, const QuantumDot &qd)
{
  os << qd.location.transpose() << " ";
  os << qd.frequency << " ";
  os << qd.damping.first << " " << qd.damping.second << " ";
  os << qd.dipole;
  return os;
}
