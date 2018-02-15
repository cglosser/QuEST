#include "math_utils.h"

int wrapmod(const int n, const int M) { return ((n % M) + M) % M; }
std::vector<double> linspace(const double low,
                             const double high,
                             const size_t n,
                             double* const step /*= nullptr*/)
{
  std::vector<double> xs(n);
  const double dx = (high - low) / (n - 1);
  if(step) *step = dx;

  for(size_t i = 0; i < n; ++i) {
    xs[i] = low + i * dx;
  }

  return xs;
}

Eigen::Vector3d unit_normal(double theta, double phi)
{
  Eigen::Vector3d rhat(std::sin(theta) * std::cos(phi),
                       std::sin(theta) * std::sin(phi), std::cos(theta));

  return rhat;
}

double gaussian(const double t) { return std::exp(-std::pow(t, 2) / 2); }
double skew_gaussian(const double alpha, const double t)
{
  return gaussian(t) * std::erfc(-alpha * t / std::sqrt(2));
}

int grid_sequence(const int n)
{
  return (1 - std::pow(-1, n) * (1 + 2 * n)) / 4;
}

std::pair<int, double> split_double(const double delay)
{
  std::pair<int, double> result;

  double idelay;
  result.second = std::modf(delay, &idelay);
  result.first = static_cast<int>(idelay);

  return result;
}

double falling_factorial(const double x, const int n)
{
  if(n == 0) return 1;
  double result = x;
  for(int i = 1; i < n; ++i) result *= x - i;

  return result;
}

std::vector<double> chebyshev_points(const int n)
{
  std::vector<double> pts(n + 1);
  for(int i = 0; i <= n; ++i) {
    pts[n - i] = (std::cos(i * M_PI / n) + 1) / 2;
  }

  return pts;
}
