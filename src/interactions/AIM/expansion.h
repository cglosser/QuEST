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

    using SpatialExpansionFunction = std::function<Eigen::Vector3cd(
        const Eigen::Vector3cd &, const Expansions::DerivativeTable &)>;

    const SpatialExpansionFunction Derivative0 = [](
        const Eigen::Vector3cd &grid_field,
        const Expansions::DerivativeTable &derivs) {
      using namespace Expansions::enums;
      return derivs[D_0] * grid_field;
    };

    const SpatialExpansionFunction GradDiv = [](
        const Eigen::Vector3cd &grid_field,
        const Expansions::DerivativeTable &derivs) {
      using namespace Expansions::enums;
      Eigen::Map<const Eigen::Matrix3d> del_del(&derivs[D_XX]);
      return del_del * grid_field;
    };

    using ExpansionFunction =
        std::function<Eigen::Vector3cd(const spacetime::vector3d<cmplx> &,
                                       const std::array<int, 4> &,
                                       const Expansions::Expansion &)>;

    //class Dyadic {
     //public:
      //Dyadic(int circulant_dimensions, double c)
          //: order(5),
            //circulant_dimensions(circulant_dimensions),
            //c(c),
            //dt_sq_coefs({3.750000000000000, -12.83333333333333,
                         //17.83333333333333, -13.00000000000000,
                         //5.083333333333333, -0.8333333333333333}){};

      //Eigen::Vector3cd operator()(const spacetime::vector3d<cmplx> &obs,
                                  //const std::array<int, 4> &coord,
                                  //const Expansions::Expansion &e)
      //{
        //Eigen::Vector3cd total_field = Eigen::Vector3cd::Zero();

        //for(int h = 0; h <= order; ++h) {
          //int w = std::max(coord[0] - h, 0) % circulant_dimensions;
          //Eigen::Map<const Eigen::Vector3cd> field(
              //&obs[w][coord[1]][coord[2]][coord[3]][0]);
          //total_field -= dt_sq_coefs[h] * Derivative0(field, e.weights);
        //}

        //total_field +=
            //std::pow(c, 2) *
            //GradDiv(Eigen::Map<const Eigen::Vector3cd>(
                        //&obs[coord[0]][coord[1]][coord[2]][coord[3]][0]),
                    //e.weights);

        //return total_field;
      //}

     //private:
      //int order, circulant_dimensions;
      //double c;
      //std::array<double, 6> dt_sq_coefs;
    //};

    auto Retardation = [](const spacetime::vector3d<cmplx> &obs,
                          const std::array<int, 4> &coord,
                          const Expansions::Expansion &e) {
      
      
    }



    class Retardation {
     public:
      Eigen::Vector3cd operator()(const spacetime::vector3d<cmplx> &obs,
                                  const std::array<int, 4> &coord,
                                  const Expansions::Expansion &e)
      {
        Eigen::Vector3cd total_field = Eigen::Vector3cd::Zero();

        for(int h = 0; h <= order; ++h) {
          int w = std::max(coord[0] - h, 0) % circulant_dimensions;
          Eigen::Map<const Eigen::Vector3cd> field(
              &obs[w][coord[1]][coord[2]][coord[3]][0]);
          total_field -= dt_coefs[h] * Derivative0(field, e.weights);
        }

        return total_field;
      }
    };

    class TimeDerivative {
     public:
      TimeDerivative(int circulant_dimensions, double c)
          : order(3),
            circulant_dimensions(circulant_dimensions),
            c(c),
            dt_coefs({1.833333333333333, -3.000000000000000, 1.500000000000000,
                      -0.3333333333333333}){};

      Eigen::Vector3cd operator()(const spacetime::vector3d<cmplx> &obs,
                                  const std::array<int, 4> &coord,
                                  const Expansions::Expansion &e)
      {
        Eigen::Vector3cd total_field = Eigen::Vector3cd::Zero();

        for(int h = 0; h <= order; ++h) {
          int w = std::max(coord[0] - h, 0) % circulant_dimensions;
          Eigen::Map<const Eigen::Vector3cd> field(
              &obs[w][coord[1]][coord[2]][coord[3]][0]);
          total_field -= dt_coefs[h] * Derivative0(field, e.weights);
        }

        return total_field;
      }

     private:
      int order, circulant_dimensions;
      double c;
      std::array<double, 6> dt_coefs;
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
