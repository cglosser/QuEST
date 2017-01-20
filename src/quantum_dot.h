#ifndef QUANTUM_DOT_H
#define QUANTUM_DOT_H

#include <Eigen/Dense>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "common.h"

class QuantumDot;

typedef std::vector<QuantumDot> DotVector;
typedef Eigen::Vector2cd matrix_elements;

class QuantumDot {
 public:
  QuantumDot() = default;
  QuantumDot(const Eigen::Vector3d &, const double,
             const std::pair<double, double> &, const Eigen::Vector3d &);

  matrix_elements liouville_rhs(const matrix_elements &, const cmplx,
                                const double) const;

  const Eigen::Vector3d &position() const { return pos; }
  const Eigen::Vector3d &dipole() const { return dip; }
  friend Eigen::Vector3d separation(const QuantumDot &, const QuantumDot &);
  friend inline double dyadic_product(const QuantumDot &obs,
                                      const Eigen::Matrix3d &dyad,
                                      const QuantumDot &src)
  {
    return obs.dip.transpose() * dyad * src.dip;
  }

  friend std::ostream &operator<<(std::ostream &, const QuantumDot &);
  friend std::istream &operator>>(std::istream &, QuantumDot &);

 private:
  Eigen::Vector3d pos;
  double freq;
  std::pair<double, double> damping;
  Eigen::Vector3d dip;
};

DotVector import_dots(const std::string &);

std::vector<
    std::function<matrix_elements(const matrix_elements &, const cmplx)>>
rhs_functions(const DotVector &, const double);

#endif
