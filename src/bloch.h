#ifndef BLOCH_H
#define BLOCH_H

#include <utility>
#include <vector>
#include <Eigen/Dense>

class QuantumDot {
 public:
  QuantumDot(const double, const double, const std::pair<double, double> &, 
      const Eigen::Vector3d &);

 private:
  double frequency, dipole;
  std::pair<double, double> damping;
  Eigen::Vector3d location;
  std::vector<Eigen::Vector3d> history; //!< Pseudospin history. The last entry represents "now" and last - i represents i steps into the past.
  Eigen::Matrix3d omega_matrix(const double, const Eigen::Vector3d &);
};

#endif
