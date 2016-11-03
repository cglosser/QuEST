#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include "lagrange_set.h"
#include "math_utils.h"

std::vector<double> input_signal();
std::vector<double> output_signal();
std::pair<int, double> delay_parts(const double);

#endif
