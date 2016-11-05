#include "bloch.h"

QuantumDot::QuantumDot():
  pos(Eigen::Vector3d(0,0,0)),
  frequency(0),
  damping(std::pair<double, double>(0,0)),
  dipole(0),
  dir(Eigen::Vector3d(0,0,0))
{}

QuantumDot::QuantumDot(
  const Eigen::Vector3d &loc,
  const double omega,
  const std::pair<double, double> &ts,
  const double dipole_strength,
  const Eigen::Vector3d &dipole_direction
):
  pos(loc),
  frequency(omega),
  damping(ts),
  dipole(dipole_strength),
  dir(dipole_direction)
{
  history.reserve(512);
}

matrix_element QuantumDot::interpolate(const UniformLagrangeSet &delay,
    int offset = 0)
{
  matrix_element result = 0;

  const size_t start = history.size() - 1 - offset;
  assert(start - config.interpolation_order >= 0);

  for(int i = 0; i <= config.interpolation_order; ++i) {
    result += history[start - i]*delay.weights[1][i];
  }

  return result;
}

std::ostream &operator<<(std::ostream &os, const QuantumDot &qd)
{
  os << qd.pos.transpose() << " ";
  os << qd.frequency << " ";
  os << qd.damping.first << " " << qd.damping.second << " ";
  os << qd.dipole << " ";
  os << qd.dir.transpose();
  return os;
}

std::istream &operator>>(std::istream &is, QuantumDot &qd)
{
  is >> qd.pos[0] >> qd.pos[1] >> qd.pos[2] >> qd.frequency
     >> qd.damping.first >> qd.damping.second >> qd.dipole
     >> qd.dir[0] >> qd.dir[1] >> qd.dir[2];

  return is;
}
