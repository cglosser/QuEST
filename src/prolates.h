#ifndef PROLATES_H
#define PROLATES_H

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <iostream>

constexpr double TOLER = 1e-5;

class Prolate {
  double alpha;
  int width;

  double sqrt_term(double);
  double d1_sqrt_term(double);
  double d2_sqrt_term(double);
  double d3_sqrt_term(double);


 public:
  Prolate(int);

  double d0(double); //0th derivative
  double d1(double); //1st derivative
  double d2(double); //2nd derivative
  double d3(double); //3rd derivative
};

#endif
