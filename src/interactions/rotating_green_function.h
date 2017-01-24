#ifndef ROTATING_GREEN_FUNCTION_H
#define ROTATING_GREEN_FUNCTION_H

#include "green_function.h"

namespace GreenFunction {
  class RotatingDyadic;
}

class GreenFunction::RotatingDyadic : public GreenFunction::Dyadic {
 public:
  RotatingDyadic(const double mu0, const double c, const double hbar, const double omega)
      : Dyadic(mu0, c, hbar), omega_(omega){};
  virtual std::vector<Eigen::Matrix3cd> coefficients(
      const Eigen::Vector3d &, const Interpolation::UniformLagrangeSet &) const;

  virtual cmplx polarization_prefactor(
      const History::soltype &matrix_elements) const
  {
    return matrix_elements[1];
  }

 private:
  double omega_;
};

#endif
