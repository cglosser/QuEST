#include "bloch.h"

QuantumDot::QuantumDot(const double omega, const double kappa, 
    const std::pair<double, double> &ts, const Eigen::Vector3d &loc) :
  history(500),
  frequency(omega),
  dipole(kappa),
  damping(ts),
  location(loc)
{}
