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
  static Array2<Expansion> get_expansions(const int,
                                          const Grid &,
                                          const std::vector<QuantumDot> &);
  Array2<Expansion> table(const std::vector<QuantumDot> &) const;
  Eigen::VectorXd q_vector(const Eigen::Vector3d &,
                           const std::array<int, 3> &) const;
  Eigen::MatrixXd w_matrix(const Eigen::Vector3d &) const;
  Eigen::VectorXd solve_expansion_system(const Eigen::Vector3d &) const;

 private:
  LeastSquaresExpansionSolver(const int box_order, const Grid &grid)
      : box_order(box_order), num_pts(std::pow(box_order + 1, 3)), grid(grid){};
  const int box_order, num_pts;
  const Grid &grid;
};

#endif
