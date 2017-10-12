#ifndef EXPANSION_H
#define EXPANSION_H

#include "grid.h"
#include "quantum_dot.h"

namespace AIM {
  struct Expansion {
    size_t index;
    double weight;
  };
  class LeastSquaresExpansionSolver;
}

class AIM::LeastSquaresExpansionSolver {
 public:
  LeastSquaresExpansionSolver(const int box_order, const Grid &grid)
      : box_order(box_order), num_pts(std::pow(box_order + 1, 3)), grid(grid){};
  Array<Expansion> table(const std::vector<QuantumDot> &) const;
  Eigen::VectorXd q_vector(const Eigen::Vector3d &) const;
  Eigen::MatrixXd w_matrix(const Eigen::Vector3d &) const;
  Eigen::VectorXd solve_expansion_system(const Eigen::Vector3d &) const;

 private:
  const int box_order, num_pts;
  const Grid &grid;
};

#endif
