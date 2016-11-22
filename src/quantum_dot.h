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

  matrix_elements liouville_rhs(const matrix_elements &, const double) const;

  friend Eigen::Vector3d separation(const QuantumDot &, const QuantumDot &);
  friend double dyadic_product(const QuantumDot &, const Eigen::Matrix3d &,
                               const QuantumDot &);

  friend std::ostream &operator<<(std::ostream &, const QuantumDot &);
  friend std::istream &operator>>(std::istream &, QuantumDot &);

 private:
  Eigen::Vector3d pos;
  double frequency;
  std::pair<double, double> damping;
  Eigen::Vector3d dipole_moment;
};

typedef std::vector<QuantumDot> DotTable;

double polarization(const matrix_elements &);

DotTable import_dots(const std::string &);

#endif
