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
  QuantumDot();
  QuantumDot(const Eigen::Vector3d &, const double,
             const std::pair<double, double> &, const Eigen::Vector3d &);


  friend std::ostream &operator<<(std::ostream &, const QuantumDot &);
  friend std::istream &operator>>(std::istream &, QuantumDot &);

  // private:
  Eigen::Vector3d pos;
  double frequency;
  std::pair<double, double> damping;
  Eigen::Vector3d dipole;
};

std::vector<QuantumDot> import_dots(const std::string &);

#endif
