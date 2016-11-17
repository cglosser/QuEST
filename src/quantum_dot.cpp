#include "quantum_dot.h"

QuantumDot::QuantumDot()
    : pos(Eigen::Vector3d(0, 0, 0)),
      frequency(0),
      damping(std::pair<double, double>(0, 0)),
      dipole(Eigen::Vector3d(0, 0, 0))
{
}

QuantumDot::QuantumDot(const Eigen::Vector3d &loc, const double omega,
                       const std::pair<double, double> &ts,
                       const Eigen::Vector3d &dip)
    : pos(loc), frequency(omega), damping(ts), dipole(dip)
{
  history.reserve(2048);
}

Eigen::Vector3d QuantumDot::polarization(const size_t idx) const
{
  const double coef = 2 * (history.at(idx)(1)).real();
  return coef * dipole;
}

std::ostream &operator<<(std::ostream &os, const QuantumDot &qd)
{
  os << qd.pos.transpose() << " ";
  os << qd.frequency << " ";
  os << qd.damping.first << " " << qd.damping.second << " ";
  os << qd.dipole.transpose();
  return os;
}

std::istream &operator>>(std::istream &is, QuantumDot &qd)
{
  is >> qd.pos[0] >> qd.pos[1] >> qd.pos[2] >>
      qd.frequency >> qd.damping.first >> qd.damping.second >>
      qd.dipole[0] >> qd.dipole[1] >> qd.dipole[2];

  return is;
}
