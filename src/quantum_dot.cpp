#include "quantum_dot.h"

QuantumDot::QuantumDot(const Eigen::Vector3d &loc, const double omega,
                       const std::pair<double, double> &ts,
                       const Eigen::Vector3d &dipole)
    : pos(loc), frequency(omega), damping(ts), dipole_moment(dipole)
{
}

matrix_elements QuantumDot::liouville_rhs(const matrix_elements &rho,
                                          const double rabi) const
{
  constexpr std::complex<double> iu(0, 1);
  double m0 = -2 * rho[1].imag() * rabi;
  std::complex<double> m1(
      iu * ((2 * rho[0].real() - 1) * rabi + rho[1] * frequency));
  return matrix_elements(m0, m1);
}

Eigen::Vector3d separation(const QuantumDot &d1, const QuantumDot &d2)
{
  return d2.pos - d1.pos;
}

std::ostream &operator<<(std::ostream &os, const QuantumDot &qd)
{
  os << qd.pos.transpose() << " " << qd.frequency << " " << qd.damping.first
     << " " << qd.damping.second << " " << qd.dipole_moment.transpose();
  return os;
}

std::istream &operator>>(std::istream &is, QuantumDot &qd)
{
  is >> qd.pos[0] >> qd.pos[1] >> qd.pos[2] >> qd.frequency >>
      qd.damping.first >> qd.damping.second >> qd.dipole_moment[0] >>
      qd.dipole_moment[1] >> qd.dipole_moment[2];
  return is;
}

std::vector<QuantumDot> import_dots(const std::string &fname)
{
  std::ifstream ifs(fname);
  std::istream_iterator<QuantumDot> in_iter(ifs), eof;
  return std::vector<QuantumDot>(in_iter, eof);
}
