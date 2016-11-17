#ifndef QUANTUM_DOT_H
#define QUANTUM_DOT_H

#include <Eigen/Dense>
#include <complex>
#include <istream>
#include <ostream>
#include <utility>
#include <vector>

typedef Eigen::Vector2cd matrix_elements;

class QuantumDot {
 public:
  QuantumDot();
  QuantumDot(const Eigen::Vector3d &, const double,
             const std::pair<double, double> &, const Eigen::Vector3d &);

  std::vector<matrix_elements> history;

  Eigen::Vector3d polarization(const size_t) const;

  friend std::ostream &operator<<(std::ostream &, const QuantumDot &);
  friend std::istream &operator>>(std::istream &, QuantumDot &);

  // private:
  Eigen::Vector3d pos;
  double frequency;
  std::pair<double, double> damping;
  Eigen::Vector3d dipole;
};

#endif
