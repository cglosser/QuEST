#ifndef MATH_UTILS_H
#define MATH_UTILS_H

#include <cmath>
#include <Eigen/Dense>
#include <utility>

typedef std::pair<int, double> QuotientRemainder;

Eigen::Vector3d unit_normal(const double, const double);
double gaussian(const double, const double, const double);
QuotientRemainder split_real(const double, const double);

#endif
