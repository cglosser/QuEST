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

    class Retardation {
     public:
      Retardation(int history_length) : history_length(history_length){};
      Eigen::Vector3cd operator()(const spacetime::vector3d<cmplx> &obs,
                                  const std::array<int, 4> &coord,
                                  const Expansions::Expansion &e)
      {
        Eigen::Map<const Eigen::Vector3cd> field(
            &obs[coord[0] % history_length][coord[1]][coord[2]][coord[3]][0]);
        return spatial::Derivative0(field, e.weights);
      }

     private:
      int history_length;
    };

    class TimeDerivative {
     public:
      TimeDerivative(const int history_length, const double dt)
          : dt_coefs({{25.0 / 12, -4.0, 3.0, -4.0 / 3, 1.0 / 4}}),
            history_length(history_length)
      {
        for(auto &coef : dt_coefs) {
          coef /= dt;
        }
      };
      Eigen::Vector3cd operator()(const spacetime::vector3d<cmplx> &obs,
                                  const std::array<int, 4> &coord,
                                  const Expansions::Expansion &e)
      {
        Eigen::Vector3cd total_field = Eigen::Vector3cd::Zero();

        for(int h = 0; h < static_cast<int>(dt_coefs.size()); ++h) {
          int w = wrap_index(std::max(coord[0] - h, 0));
          Eigen::Map<const Eigen::Vector3cd> field(
              &obs[w][coord[1]][coord[2]][coord[3]][0]);
          total_field += dt_coefs[h] * spatial::Derivative0(field, e.weights);
        }
        return total_field;
      }

     private:
      std::array<double, 5> dt_coefs;
      int history_length;
      int wrap_index(int t) { return t % history_length; };
    };

    class EFIE {
     public:
      EFIE(int history_length, double c, double dt)
          : dt2_coefs(
                {{15.0 / 4, -77.0 / 6, 107.0 / 6, -13.0, 61.0 / 12, -5.0 / 6}}),
            history_length(history_length),
            c(c)
      {
        for(auto &coef : dt2_coefs) {
          coef /= std::pow(dt, 2);
        }
      };

      Eigen::Vector3cd operator()(const spacetime::vector3d<cmplx> &obs,
                                  const std::array<int, 4> &coord,
                                  const Expansions::Expansion &e)
      {
        Eigen::Vector3cd dt2 = Eigen::Vector3cd::Zero();
        for(int h = 0; h < static_cast<int>(dt2_coefs.size()); ++h) {
          int w = wrap_index(std::max(coord[0] - h, 0));
          Eigen::Map<const Eigen::Vector3cd> field(
              &obs[w][coord[1]][coord[2]][coord[3]][0]);
          dt2 += dt2_coefs[h] * spatial::Derivative0(field, e.weights);
        }

        Eigen::Map<const Eigen::Vector3cd> field(
            &obs[wrap_index(coord[0])][coord[1]][coord[2]][coord[3]][0]);
        Eigen::Vector3cd del_del = spatial::GradDiv(field, e.weights);

        return dt2 - std::pow(c, 2) * del_del;
      }

     private:
      std::array<double, 6> dt2_coefs;
      int history_length;
      double c;
      int wrap_index(int t) { return t % history_length; };
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
