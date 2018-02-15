#ifndef EXPANSION_H
#define EXPANSION_H

#include "grid.h"
#include "quantum_dot.h"
#include "spacetime.h"

namespace AIM {
  namespace Expansions {
    struct Expansion {
      size_t index;
      double weight;
    };

    using ExpansionTable = boost::multi_array<Expansion, 2>;
    class LeastSquaresExpansionSolver;
  }
}

class AIM::Expansions::LeastSquaresExpansionSolver {
 public:
  static ExpansionTable get_expansions(const int,
                                       const Grid &,
                                       const std::vector<QuantumDot> &);
  Eigen::VectorXd solve(const Eigen::Vector3d &) const;
  ExpansionTable table(const std::vector<QuantumDot> &) const;
  Eigen::MatrixXd w_matrix(const Eigen::Vector3d &) const;

 private:
  int box_order, num_pts;
  Eigen::VectorXd qvec;
  const Grid &grid;

  LeastSquaresExpansionSolver(const int box_order, const Grid &grid)
      : box_order(box_order),
        num_pts(std::pow(box_order + 1, 3)),
        qvec(Eigen::VectorXd::Zero(num_pts)),
        grid(grid)
  {
    qvec(0) = 1;
  };
};

#endif
