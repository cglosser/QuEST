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
  QuantumDot(const Eigen::Vector3d &, const double,
      const std::pair<double, double> &, const double, const Eigen::Vector3d &);

  std::vector<matrix_element> history;

  friend std::ostream &operator<<(std::ostream &, const QuantumDot &);

 //private:
  Eigen::Vector3d location;
  double frequency;
  std::pair<double, double> damping;
  double dipole;
  Eigen::Vector3d orientation;
};

#endif
