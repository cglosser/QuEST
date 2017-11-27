#ifndef EXPANSION_H
#define EXPANSION_H

#include "grid.h"
#include "quantum_dot.h"

namespace AIM {
  namespace Expansions {
    inline namespace enums {
      constexpr int NUM_DERIVS = 13;
      enum DERIVATIVE_ORDER {
        D_0,
        D_X,
        D_Y,
        D_Z,
        D_XX,
        D_XY,
        D_XZ,
        D_YX,
        D_YY,
        D_YZ,
        D_ZX,
        D_ZY,
        D_ZZ
      };
    }

    using DerivativeTable = std::array<double, enums::NUM_DERIVS>;

    struct Expansion {
      size_t index;
      DerivativeTable weights;
    };
    using ExpansionTable = boost::multi_array<Expansion, 2>;
    class LeastSquaresExpansionSolver;

    using ExpansionFunction = std::function<Eigen::Vector3cd(
        const Eigen::Vector3cd &, const Expansions::DerivativeTable &)>;

    const ExpansionFunction Identity = [](
        const Eigen::Vector3cd &grid_field,
        const Expansions::DerivativeTable &derivs) {
      using namespace Expansions::enums;
      return derivs[D_0] * grid_field;
    };

    const ExpansionFunction DerivX = [](
        const Eigen::Vector3cd &grid_field,
        const Expansions::DerivativeTable &derivs) {
      using namespace Expansions::enums;
      return derivs[D_X] * grid_field;
    };

    const ExpansionFunction DerivXX = [](
        const Eigen::Vector3cd &grid_field,
        const Expansions::DerivativeTable &derivs) {
      using namespace Expansions::enums;
      return derivs[D_XX] * grid_field;
    };

    const ExpansionFunction GradDiv = [](
        const Eigen::Vector3cd &grid_field,
        const Expansions::DerivativeTable &derivs) {
      using namespace Expansions::enums;
      Eigen::Map<const Eigen::Matrix3d> del_del(&derivs[D_XX]);
      return del_del * grid_field;
      // return grid_field;
    };
  }
}

class AIM::Expansions::LeastSquaresExpansionSolver {
 public:
  static ExpansionTable get_expansions(const int,
                                       const Grid &,
                                       const std::vector<QuantumDot> &);
  ExpansionTable table(const std::vector<QuantumDot> &) const;
  Eigen::VectorXd q_vector(const std::array<int, 3> & = {{0, 0, 0}}) const;
  Eigen::MatrixXd w_matrix(const Eigen::Vector3d &) const;

 private:
  LeastSquaresExpansionSolver(const int box_order, const Grid &grid)
      : box_order(box_order), num_pts(std::pow(box_order + 1, 3)), grid(grid){};
  const int box_order, num_pts;
  const Grid &grid;
};

#endif
