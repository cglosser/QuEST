#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <Eigen/Dense>
#include <cmath>
#include <utility>
#include <vector>

std::vector<double> linspace(const double,
                             const double,
                             const size_t,
                             double* const = nullptr);
Eigen::Vector3d unit_normal(const double, const double);
double gaussian(const double);
double skew_gaussian(const double, const double);
int grid_sequence(const size_t);
std::pair<int, double> split_double(const double);
double falling_factorial(const double, int);
#endif
