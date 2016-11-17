#ifndef QUANTUM_DOT_H
#define QUANTUM_DOT_H

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <complex>
#include <istream>
#include <ostream>
#include <utility>
#include <vector>

typedef Eigen::Vector2cd matrix_elements;
typedef std::vector<matrix_elements, Eigen::aligned_allocator<matrix_elements>>
    history_vector;

class QuantumDot {
 public:
  QuantumDot();
  QuantumDot(const Eigen::Vector3d &, const double,
             const std::pair<double, double> &, const Eigen::Vector3d &);

  history_vector history;

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
