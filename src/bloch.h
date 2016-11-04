#ifndef BLOCH_H
#define BLOCH_H

#include <complex>
#include <ostream>
#include <utility>
#include <vector>
#include <Eigen/Dense>

typedef std::pair<double, std::complex<double>> matrix_element;

class QuantumDot {
 public:
  QuantumDot(const double, const double, const std::pair<double, double> &, 
      const Eigen::Vector3d &);

  std::vector<matrix_element> history;

  friend std::ostream &operator<<(std::ostream &, const QuantumDot &);

 private:
  double frequency, dipole;
  std::pair<double, double> damping;
  Eigen::Vector3d location;
};

#endif
