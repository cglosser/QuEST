#ifndef QUANTUM_DOT_H
#define QUANTUM_DOT_H

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <complex>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>

typedef Eigen::Vector2cd matrix_elements;

class QuantumDot {
 public:
  QuantumDot() = default;
  QuantumDot(const Eigen::Vector3d &, const double,
             const std::pair<double, double> &, const Eigen::Vector3d &);

  matrix_elements liouville_rhs(const matrix_elements &, const double) const;

  const Eigen::Vector3d &position() const { return pos; }
  const Eigen::Vector3d &dipole() const { return dipole_moment; }

  friend Eigen::Vector3d separation(const QuantumDot &, const QuantumDot &);
  friend inline double dyadic_product(const QuantumDot &obs,
                                      const Eigen::Matrix3d &dyad,
                                      const QuantumDot &src)
  {
    return obs.dipole_moment.transpose() * dyad * src.dipole_moment;
  }

  friend std::ostream &operator<<(std::ostream &, const QuantumDot &);
  friend std::istream &operator>>(std::istream &, QuantumDot &);

 private:
  Eigen::Vector3d pos;
  double frequency;
  std::pair<double, double> damping;
  Eigen::Vector3d dipole_moment;
};

std::vector<QuantumDot> import_dots(const std::string &);

inline double polarization(const matrix_elements &mel)
{
  return 2 * mel[1].real();
}

#endif
