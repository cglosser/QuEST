#include "math_utils.h"

int Math::wrapmod(const int n, const int M) { return ((n % M) + M) % M; }
std::vector<double> Math::linspace(const double low,
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

Eigen::Vector3d Math::unit_normal(double theta, double phi)
{
  Eigen::Vector3d rhat(std::sin(theta) * std::cos(phi),
                       std::sin(theta) * std::sin(phi), std::cos(theta));

  return rhat;
}

double Math::gaussian(const double t) { return std::exp(-std::pow(t, 2) / 2); }
double Math::skew_gaussian(const double alpha, const double t)
{
  return gaussian(t) * std::erfc(-alpha * t / std::sqrt(2));
}

int Math::grid_sequence(const int n)
{
  return (1 - std::pow(-1, n) * (1 + 2 * n)) / 4;
}

std::pair<int, double> Math::split_double(const double delay)
{
  std::pair<int, double> result;

  double idelay;
  result.second = std::modf(delay, &idelay);
  result.first = static_cast<int>(idelay);

  return result;
}

double Math::falling_factorial(const double x, const int n)
{
  if(n == 0) return 1;
  double result = x;
  for(int i = 1; i < n; ++i) result *= x - i;

  return result;
}

std::vector<double> Math::chebyshev_points(const int n)
{
  std::vector<double> pts(n + 1);
  for(int i = 0; i <= n; ++i) {
    pts[n - i] = (std::cos(i * M_PI / n) + 1) / 2;
  }

  return pts;
}

double Math::Chebyshev::T(int n, double x)
{
  switch(n) {
    case(0): return T0(x);
    case(1): return T1(x);
    case(2): return T2(x);
    case(3): return T3(x);
    case(4): return T4(x);
    case(5): return T5(x);
    case(6): return T6(x);
    case(7): return T7(x);
    case(8): return T8(x);
  }

  return 0;
}
constexpr double Math::Chebyshev::T0(__attribute__((unused)) const double x)
{
  return 1;
}
constexpr double Math::Chebyshev::T1(const double x) { return x; }
constexpr double Math::Chebyshev::T2(const double x) { return -1 + 2 * x * x; }
constexpr double Math::Chebyshev::T3(const double x)
{
  return x * (-3 + 4 * x * x);
}
constexpr double Math::Chebyshev::T4(const double x)
{
  return 1 + x * x * (-8 + 8 * x * x);
}
constexpr double Math::Chebyshev::T5(const double x)
{
  return x * (5 + x * x * (-20 + 16 * x * x));
}
constexpr double Math::Chebyshev::T6(const double x)
{
  return -1 + x * x * (18 + x * x * (-48 + 32 * x * x));
}
constexpr double Math::Chebyshev::T7(const double x)
{
  return x * (-7 + x * x * (56 + x * x * (-112 + 64 * x * x)));
}
constexpr double Math::Chebyshev::T8(const double x)
{
  return 1 + x * x * (-32 + x * x * (160 + x * x * (-256 + 128 * x * x)));
}

double Math::Chebyshev::T_d1(int n, double x)
{
  switch(n) {
    case(0): return T0_d1(x);
    case(1): return T1_d1(x);
    case(2): return T2_d1(x);
    case(3): return T3_d1(x);
    case(4): return T4_d1(x);
    case(5): return T5_d1(x);
    case(6): return T6_d1(x);
    case(7): return T7_d1(x);
    case(8): return T8_d1(x);
  }

  return 0;
}
constexpr double Math::Chebyshev::T0_d1(__attribute__((unused)) const double x)
{
  return 0;
}
constexpr double Math::Chebyshev::T1_d1(__attribute__((unused)) const double x)
{
  return 1;
}
constexpr double Math::Chebyshev::T2_d1(const double x) { return 4 * x; }
constexpr double Math::Chebyshev::T3_d1(const double x)
{
  return -3 + 12 * x * x;
}
constexpr double Math::Chebyshev::T4_d1(const double x)
{
  return x * (-16 + 32 * x * x);
}
constexpr double Math::Chebyshev::T5_d1(const double x)
{
  return 5 + x * x * (-60 + 80 * x * x);
}
constexpr double Math::Chebyshev::T6_d1(const double x)
{
  return x * (36 + x * x * (-192 + 192 * x * x));
}
constexpr double Math::Chebyshev::T7_d1(const double x)
{
  return -7 + x * x * (168 + x * x * (-560 + 448 * x * x));
}
constexpr double Math::Chebyshev::T8_d1(const double x)
{
  return x * (-64 + x * x * (640 + x * x * (-1536 + 1024 * x * x)));
}

double Math::Chebyshev::T_d2(int n, double x)
{
  switch(n) {
    case(0): return T0_d2(x);
    case(1): return T1_d2(x);
    case(2): return T2_d2(x);
    case(3): return T3_d2(x);
    case(4): return T4_d2(x);
    case(5): return T5_d2(x);
    case(6): return T6_d2(x);
    case(7): return T7_d2(x);
    case(8): return T8_d2(x);
  }

  return 0;
}
constexpr double Math::Chebyshev::T0_d2(__attribute__((unused)) const double x)
{
  return 0;
}
constexpr double Math::Chebyshev::T1_d2(__attribute__((unused)) const double x)
{
  return 0;
}
constexpr double Math::Chebyshev::T2_d2(__attribute__((unused)) const double x)
{
  return 4;
}
constexpr double Math::Chebyshev::T3_d2(const double x) { return 24 * x; }
constexpr double Math::Chebyshev::T4_d2(const double x)
{
  return -16 + 96 * x * x;
}
constexpr double Math::Chebyshev::T5_d2(const double x)
{
  return x * (-120 + 320 * x * x);
}
constexpr double Math::Chebyshev::T6_d2(const double x)
{
  return 36 + x * x * (-576 + 960 * x * x);
}
constexpr double Math::Chebyshev::T7_d2(const double x)
{
  return x * (336 + x * x * (-2240 + 26882 * x * x));
}
constexpr double Math::Chebyshev::T8_d2(const double x)
{
  return -64 + x * x * (1920 + x * x * (-7680 + 7168 * x * x));
}
