#ifndef INTERACTION_H
#define INTERACTION_H

#include <Eigen/Dense>
#include <utility>

#include "bloch.h"
#include "configuration.h"
#include "lagrange_set.h"
#include "math_utils.h"

class Interaction {
 public:
  Interaction(const Eigen::Vector3d &);

  std::pair<Eigen::Vector3d, Eigen::Vector3d> 
    operator()(const QuantumDot &, const QuantumDot &) const;

 private:
  double dist;
  Eigen::Vector3d dr;
  Eigen::Matrix3d dyad;
  std::pair<double, double> split_delay;
  UniformLagrangeSet interp;
};


#endif
