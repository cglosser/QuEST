#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <Eigen/Dense>
#include <cmath>
#include <utility>
#include <vector>

namespace Math {
  int wrapmod(const int, const int);
  std::vector<double> linspace(const double,
                               const double,
                               const size_t,
                               double* const = nullptr);
  Eigen::Vector3d unit_normal(const double, const double);
  double gaussian(const double);
  double skew_gaussian(const double, const double);
  int grid_sequence(const int);
  std::pair<int, double> split_double(const double);
  double falling_factorial(const double, int);
  std::vector<double> chebyshev_points(const int);
  double ChebyshevT(const int, const double);

  namespace Chebyshev {
    double T(int, double);
    constexpr double T0(const double);
    constexpr double T1(const double);
    constexpr double T2(const double);
    constexpr double T3(const double);
    constexpr double T4(const double);
    constexpr double T5(const double);
    constexpr double T6(const double);
    constexpr double T7(const double);
    constexpr double T8(const double);

    double T_d1(int, double);
    constexpr double T0_d1(const double);
    constexpr double T1_d1(const double);
    constexpr double T2_d1(const double);
    constexpr double T3_d1(const double);
    constexpr double T4_d1(const double);
    constexpr double T5_d1(const double);
    constexpr double T6_d1(const double);
    constexpr double T7_d1(const double);
    constexpr double T8_d1(const double);

    double T_d2(int, double);
    constexpr double T0_d2(const double);
    constexpr double T1_d2(const double);
    constexpr double T2_d2(const double);
    constexpr double T3_d2(const double);
    constexpr double T4_d2(const double);
    constexpr double T5_d2(const double);
    constexpr double T6_d2(const double);
    constexpr double T7_d2(const double);
    constexpr double T8_d2(const double);
  }
}

#endif
