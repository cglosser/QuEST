#include "math_utils.h"

Eigen::Vector3d unit_normal(double theta, double phi) {
  Eigen::Vector3d rhat(
    std::sin(theta)*std::cos(phi),
    std::sin(theta)*std::sin(phi),
    std::cos(theta)
  );

  return rhat;
}

double gaussian(const double t, const double mu, const double sigma)
{
  return std::exp(-std::pow((t - mu)/sigma, 2)/2);
}

QuotientRemainder split_real(const double num, const double denom)
{
  int quotient = std::floor(num/denom);
  double remainder = std::fmod(num, denom);

  return QuotientRemainder(quotient, remainder);
}
