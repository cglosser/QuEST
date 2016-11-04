#include "bloch.h"

QuantumDot::QuantumDot(
  const Eigen::Vector3d &loc,
  const double omega,
  const std::pair<double, double> &ts,
  const double dipole_strength,
  const Eigen::Vector3d &dipole_direction
):
  history(500),
  location(loc),
  frequency(omega),
  damping(ts),
  dipole(dipole_strength),
  orientation(dipole_direction)
{}

std::ostream &operator<<(std::ostream &os, const QuantumDot &qd)
{
  os << qd.location.transpose() << " ";
  os << qd.frequency << " ";
  os << qd.damping.first << " " << qd.damping.second << " ";
  os << qd.dipole << " ";
  os << qd.orientation.transpose();
  return os;
}
