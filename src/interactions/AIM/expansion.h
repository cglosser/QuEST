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

    struct Expansion {
      size_t index;
      double d0;
      Eigen::Vector3d del;
      Eigen::Matrix3d del_sq;
    };

    using ExpansionTable = boost::multi_array<Expansion, 2>;
    class LeastSquaresExpansionSolver;

    using ExpansionFunction =
        std::function<Eigen::Vector3cd(const spacetime::vector3d<cmplx> &,
                                       const std::array<int, 4> &,
                                       const Expansions::Expansion &)>;

    class RetardationBase {
     public:
      RetardationBase(const int history_length)
          : history_length(history_length){};

     protected:
      int history_length;
      int wrap_index(int t) const { return t % history_length; };
    };

    class Retardation : public RetardationBase {
     public:
      Retardation(int history_length) : RetardationBase(history_length){};
      Eigen::Vector3cd operator()(const spacetime::vector3d<cmplx> &obs,
                                  const std::array<int, 4> &coord,
                                  const Expansions::Expansion &e)
      {
        Eigen::Map<const Eigen::Vector3cd> field(
            &obs[wrap_index(coord[0])][coord[1]][coord[2]][coord[3]][0]);
        return e.d0 * field;
      }
    };

    class TimeDerivative : public RetardationBase {
     public:
      TimeDerivative(const int history_length, const double dt)
          : RetardationBase(history_length),
            dt_coefs({{25.0 / 12, -4.0, 3.0, -4.0 / 3, 1.0 / 4}})
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
          total_field += dt_coefs[h] * e.d0 * field;
        }
        return total_field;
      }

     private:
      std::array<double, 5> dt_coefs;
    };

    class Oper : public RetardationBase {
     public:
      Oper(int history_length) : RetardationBase(history_length){};

      Eigen::Vector3cd operator()(const spacetime::vector3d<cmplx> &obs,
                                  const std::array<int, 4> &coord,
                                  const Expansions::Expansion &e)
      {
        Eigen::Map<const Eigen::Vector3cd> field(
            &obs[wrap_index(coord[0])][coord[1]][coord[2]][coord[3]][0]);
        return e.del_sq * field;
      }
    };

    class Del_Del : public RetardationBase {
     public:
      Del_Del(int history_length) : RetardationBase(history_length){};

      Eigen::Vector3cd operator()(const spacetime::vector3d<cmplx> &obs,
                                  const std::array<int, 4> &coord,
                                  const Expansions::Expansion &e)
      {
        Eigen::Map<const Eigen::Vector3cd> field(
            &obs[wrap_index(coord[0])][coord[1]][coord[2]][coord[3]][0]);
        return e.del_sq * field;
      }
    };

    class EFIE : public RetardationBase {
     public:
      EFIE(int history_length, double c, double dt)
          : RetardationBase(history_length),
            dt2_coefs(
                {{15.0 / 4, -77.0 / 6, 107.0 / 6, -13.0, 61.0 / 12, -5.0 / 6}}),
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
          dt2 += dt2_coefs[h] * e.d0 * field;
        }

        Eigen::Map<const Eigen::Vector3cd> field(
            &obs[wrap_index(coord[0])][coord[1]][coord[2]][coord[3]][0]);
        Eigen::Vector3cd del_sq = e.del_sq * field;

        return -dt2 + std::pow(c, 2) * del_sq;
      }

     private:
      std::array<double, 6> dt2_coefs;
      double c;
    };

    class RotatingEFIE : public RetardationBase {
     public:
      RotatingEFIE(int history_length, double c, double dt, double omega)
          : RetardationBase(history_length),
            dt1_coefs({{137.0 / 60, -5.0, 5.0, -10.0 / 3, 5.0 / 4, -1.0 / 5}}),
            dt2_coefs(
                {{15.0 / 4, -77.0 / 6, 107.0 / 6, -13.0, 61.0 / 12, -5.0 / 6}}),
            c(c),
            omega(omega)
      {
        for(auto &coef : dt1_coefs) coef /= std::pow(dt, 1);
        for(auto &coef : dt2_coefs) coef /= std::pow(dt, 2);
      }

      Eigen::Vector3cd operator()(const spacetime::vector3d<cmplx> &obs,
                                  const std::array<int, 4> &coord,
                                  const Expansions::Expansion &e)
      {
        Eigen::Matrix3cd deriv = Eigen::Matrix3cd::Zero();
        Eigen::Map<const Eigen::Vector3cd> present_field(
            &obs[wrap_index(coord[0])][coord[1]][coord[2]][coord[3]][0]);

        deriv.col(0) = e.d0 * present_field;

        for(int h = 0; h < static_cast<int>(dt2_coefs.size()); ++h) {
          const int w = wrap_index(std::max(coord[0] - h, 0));
          Eigen::Map<const Eigen::Vector3cd> past_field(
              &obs[w][coord[1]][coord[2]][coord[3]][0]);
          deriv.col(1) += dt1_coefs[h] * e.d0 * past_field;
          deriv.col(2) += dt2_coefs[h] * e.d0 * past_field;
        }

        // Notice that this is just the derivative of E(t)exp(iwt), just
        // without the exp(iwt) because it gets suppressed in doing RWA
        // calculations. There is no exp(-ikr) factor here as that gets
        // taken care of in the normalization class.

        const auto time = deriv.col(2) + 2.0 * iu * omega * deriv.col(1) -
                          std::pow(omega, 2) * deriv.col(0);
        const auto space = e.del_sq * present_field;

        return -time + std::pow(c, 2) * space;
      }

     private:
      std::array<double, 6> dt1_coefs, dt2_coefs;
      double c, omega;
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
