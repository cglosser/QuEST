#include "quantum_dot.h"

QuantumDot::QuantumDot()
    : pos(Eigen::Vector3d(0, 0, 0)),
      frequency(0),
      damping(std::pair<double, double>(0, 0)),
      dipole_moment(Eigen::Vector3d(0, 0, 0))
{
}

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

double dyadic_product(const QuantumDot &obs, const Eigen::Matrix3d &dyad,
                      const QuantumDot &src)
{
  return obs.dipole_moment.transpose() * dyad * src.dipole_moment;
}

std::ostream &operator<<(std::ostream &os, const QuantumDot &qd)
{
  os << qd.pos.transpose() << " ";
  os << qd.frequency << " ";
  os << qd.damping.first << " " << qd.damping.second << " ";
  os << qd.dipole_moment.transpose();
  return os;
}

std::istream &operator>>(std::istream &is, QuantumDot &qd)
{
  is >> qd.pos[0] >> qd.pos[1] >> qd.pos[2] >> qd.frequency >>
      qd.damping.first >> qd.damping.second >> qd.dipole_moment[0] >>
      qd.dipole_moment[1] >> qd.dipole_moment[2];

  return is;
}

double polarization(const matrix_elements &mel) { return 2 * mel[1].real(); }
std::vector<QuantumDot> import_dots(const std::string &fname)
{
  std::ifstream ifs(fname);
  std::istream_iterator<QuantumDot> in_iter(ifs), eof;
  return std::vector<QuantumDot>(in_iter, eof);
}
