#include "math_utils.h"

std::vector<double> linspace(const double low, const double high, const size_t n,
    double * const step /*= nullptr*/)
{
  std::vector<double> xs(n);
  const double dx = (high - low)/(n - 1);
  if(step) *step = dx;

  for(size_t i = 0; i < n; ++i) {
    xs[i] = low + i*dx;
  }

  return xs;
}

Eigen::Vector3d unit_normal(double theta, double phi) {
  Eigen::Vector3d rhat(
    std::sin(theta)*std::cos(phi),
    std::sin(theta)*std::sin(phi),
    std::cos(theta)
  );

  return rhat;
}

double gaussian(const double t)
{
  return std::exp(-std::pow(t, 2)/2);
}

double skew_gaussian(const double alpha, const double t)
{
  return gaussian(t)*std::erfc(-alpha*t/std::sqrt(2));
}
