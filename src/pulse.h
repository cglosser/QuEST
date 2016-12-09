#ifndef PULSE_H
#define PULSE_H

#include <Eigen/Dense>
#include <iostream>

#include "math_utils.h"

class Pulse {
 public:
  Pulse() = default;
  Pulse(const double, const double, const double, const double,
        const Eigen::Vector3d &, const Eigen::Vector3d &);

  Eigen::Vector3d operator()(const Eigen::Vector3d &, const double);

  friend std::ostream &operator<<(std::ostream &, const Pulse &);
  friend std::istream &operator>>(std::istream &, Pulse &);

 private:
  double amplitude, delay, width, frequency;
  Eigen::Vector3d wavevector, polarization;
};

#endif
