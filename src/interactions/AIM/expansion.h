#ifndef EXPANSION_H
#define EXPANSION_H

#include "grid.h"
#include "quantum_dot.h"
#include "spacetime.h"

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

    namespace spatial {
      using ExpansionFunction = std::function<Eigen::Vector3cd(
          const Eigen::Vector3cd &, const Expansions::DerivativeTable &)>;

      const ExpansionFunction Derivative0 = [](
          const Eigen::Vector3cd &grid_field,
          const Expansions::DerivativeTable &derivs) -> Eigen::Vector3cd {
        using namespace Expansions::enums;
        return derivs[D_0] * grid_field;
      };

      const ExpansionFunction GradDiv = [](
          const Eigen::Vector3cd &grid_field,
          const Expansions::DerivativeTable &derivs) -> Eigen::Vector3cd {
        using namespace Expansions::enums;
        Eigen::Map<const Eigen::Matrix3d> del_del(&derivs[D_XX]);
        return del_del * grid_field;
      };
    }

    using ExpansionFunction =
        std::function<Eigen::Vector3cd(const spacetime::vector3d<cmplx> &,
                                       const std::array<int, 4> &,
                                       const Expansions::Expansion &)>;

    const ExpansionFunction Retardation = [](
        const spacetime::vector3d<cmplx> &obs,
        const std::array<int, 4> &coord,
        const Expansions::Expansion &e) -> Eigen::Vector3cd {
      Eigen::Map<const Eigen::Vector3cd> field(
          &obs[coord[0]][coord[1]][coord[2]][coord[3]][0]);
      return spatial::Derivative0(field, e.weights);
    };

    class TimeDerivative {
     public:
      TimeDerivative(int history_length)
          : history_length(history_length),
            dt_coefs({{1.0, 0.0, 0.0, 0.0, 0.0, 0.0}}){};
      // dt_coefs({{137./60, -5., 5., -10./3, 5./4, -1./5}}){};
      // dt_coefs({{1.833333333333333, -3.000000000000000, 1.500000000000000,
      //-0.3333333333333333}}){};

      Eigen::Vector3cd operator()(const spacetime::vector3d<cmplx> &obs,
                                  const std::array<int, 4> &coord,
                                  const Expansions::Expansion &e)
      {
        Eigen::Vector3cd total_field = Eigen::Vector3cd::Zero();

        // for(int h = 0; h < static_cast<int>(dt_coefs.size()); ++h) {
        for(int h = 0; h < 1; ++h) {
          int w = wrap_index(std::max(coord[0] - h, 0));
          Eigen::Map<const Eigen::Vector3cd> field(
              &obs[w][coord[1]][coord[2]][coord[3]][0]);
          total_field += dt_coefs[h] * spatial::Derivative0(field, e.weights);
        }

        return total_field;
      }

     private:
      int history_length;
      std::array<double, 6> dt_coefs;
      int wrap_index(int t)
      {
        return t % (history_length + 3); // FIX ME! THIS IS IN PROGRESS!
        // return (t % history_length + history_length) % history_length;
      };
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
