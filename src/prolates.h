#ifndef PROLATES_H
#define PROLATES_H

#include <cmath>

constexpr double TOLER = 1e-5;

class Prolate {
  double omega0, dt;
 public:
  Prolate(double, double);
};

#endif
