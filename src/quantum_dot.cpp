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
}

matrix_elements liouville_rhs(const matrix_elements &rho, const double rabi) const
{
  double m0 = -2 * rho[1].imag() * rabi;
  std::complex<double> m1(0, (2*rho[0] - 1)*rabi + rho[1]*frequency)
  return matrix_elements(m0, m1);
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

std::vector<QuantumDot> import_dots(const std::string &fname)
{
  std::ifstream ifs(fname);
  std::istream_iterator<QuantumDot> in_iter(ifs), eof;
  return std::vector<QuantumDot>(in_iter, eof);
}
