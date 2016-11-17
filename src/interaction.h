#ifndef INTERACTION_H
#define INTERACTION_H

#include <Eigen/Dense>
#include <utility>

#include "configuration.h"
#include "lagrange_set.h"
#include "math_utils.h"
#include "quantum_dot.h"

class Interaction {
 public:
  Interaction(const QuantumDot &, const QuantumDot &);

 //private:
  std::pair<int, double> delay;
  std::pair<int, double> compute_delay(const double) const;

  Eigen::Matrix3d rhat_dyadic(const Eigen::Vector3d &) const;
  //Eigen::Vector3d nearfield_dyadic(const Eigen::Vector3d &) const;
  //Eigen::Vector3d midfield_dyadic(const Eigen::Vector3d &) const;
  //Eigen::Vector3d farfield_dyadic(const Eigen::Vector3d &) const;
};


#endif
