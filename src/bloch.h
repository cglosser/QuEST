#ifndef BLOCH_H
#define BLOCH_H

#include <complex>
#include <istream>
#include <ostream>
#include <utility>
#include <vector>
#include <Eigen/Dense>

#include "lagrange_set.h"

//typedef Eigen::Vector2cd matrix_element;
typedef double matrix_element;

class QuantumDot {
  public:
  QuantumDot();
  QuantumDot(const Eigen::Vector3d &, const double,
      const std::pair<double, double> &, const double, const Eigen::Vector3d &);

  std::vector<matrix_element> history;
  matrix_element interpolate(const UniformLagrangeSet &, int);

  friend std::ostream &operator<<(std::ostream &, const QuantumDot &);
  friend std::istream &operator>>(std::istream &, QuantumDot &);

 //private:
  Eigen::Vector3d pos;
  double frequency;
  std::pair<double, double> damping;
  double dipole;
  Eigen::Vector3d dir;
};

#endif
