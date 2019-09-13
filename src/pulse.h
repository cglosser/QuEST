#ifndef PULSE_H
#define PULSE_H

#include <Eigen/Dense>
#include <fstream>
#include <iostream>

#include "math_utils.h"
#include "configuration.h"

class Pulse {
 public:
  Pulse() = default;
  Pulse(const double, const double, const double, const double,
        const Eigen::Vector3d &, const Eigen::Vector3d &, const Configuration::REFERENCE_FRAME); //add a fixed vs rot argument

  Eigen::Vector3d operator()(const Eigen::Vector3d &, const double) const;

  friend std::ostream &operator<<(std::ostream &, const Pulse &);
  friend std::istream &operator>>(std::istream &, Pulse &);
  friend void set_reference_frame(Pulse &, Configuration::REFERENCE_FRAME);

 private:
  double amplitude, delay, width, freq;
  Eigen::Vector3d wavevector, polarization;

  Configuration::REFERENCE_FRAME ref_frame;//add a fixed vs rot attribute
};

Pulse read_pulse_config(const std::string &, const Configuration::REFERENCE_FRAME); //add a config.ref_frame argument

#endif
