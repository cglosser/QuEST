#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <cmath>
#include <Eigen/Dense>
#include <utility>
#include <vector>

std::vector<double> linspace(const double, const double,
    const size_t, double * const = nullptr);
Eigen::Vector3d unit_normal(const double, const double);
double gaussian(const double);
double skew_gaussian(const double, const double);

#endif
