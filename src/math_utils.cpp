#include "math_utils.h"

Eigen::Vector3d unit_normal(double theta, double phi) {
  Eigen::Vector3d rhat(
    std::sin(theta)*std::cos(phi),
    std::sin(theta)*std::sin(phi),
    std::cos(theta)
  );

  return rhat;
}
