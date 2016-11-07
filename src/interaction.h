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
  Interaction(const Eigen::Vector3d &);

  std::pair<Eigen::Vector3d, Eigen::Vector3d>
    operator()(const QuantumDot &, const QuantumDot &) const;

 //private:
  double dist;
  Eigen::Vector3d dr;
  std::pair<int, double> delay;
  UniformLagrangeSet interp;

  Eigen::Matrix3d rhat_dyadic() const;
};


#endif
