#include "bloch.h"

QuantumDot::QuantumDot(const double omega, const double kappa, 
    const std::pair<double, double> &ts, const Eigen::Vector3d &loc)
{
  frequency = omega;
  dipole    = kappa;
  damping   = ts;
  location  = loc;
}

Eigen::Matrix3d QuantumDot::omega_matrix(const double field_freq, 
    const Eigen::Vector3d &efield)
{
  double detuning = frequency - field_freq;
  double rabi_freq = dipole*efield.norm();

  Eigen::Matrix3d omega;
  omega << -1/damping.second, -detuning,         0,
           detuning,          -1/damping.second, rabi_freq,
           0,                 -rabi_freq,        -1/damping.first;

  return omega;
}
