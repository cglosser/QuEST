#ifndef EXPANSION_H
#define EXPANSION_H

#include "grid.h"
#include "quantum_dot.h"
#include "spacetime.h"

namespace AIM {
  namespace Expansions {
    struct Expansion {
      int index;
      double weight;
    };

    using ExpansionTable = boost::multi_array<Expansion, 2>;
    class LeastSquaresExpansionSolver;
  }
}

class AIM::Expansions::LeastSquaresExpansionSolver {
 public:
  LeastSquaresExpansionSolver(const Grid &grid)
      : num_pts(std::pow(grid.order() + 1, 3)),
        qvec(Eigen::VectorXd::Zero(num_pts)),
        grid(grid)
  {
    qvec(0) = 1;
  };

  Eigen::MatrixXd w_matrix(const Eigen::Vector3d &) const;
  Eigen::VectorXd solve(const Eigen::Vector3d &) const;
  ExpansionTable table(const std::vector<QuantumDot> &) const;
  boost::multi_array<double, 4> chebyshev_lambda_weights(
      const std::vector<double> &) const;

 private:
  int num_pts;
  Eigen::VectorXd qvec;
  const Grid &grid;
};

#endif
