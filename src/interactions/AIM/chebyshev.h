#ifndef CHEBYSHEV_H
#define CHEBYSHEV_H

#include <array>
#include <boost/multi_array.hpp>
#include <cmath>
#include "grid.h"

template <int M>
class Chebyshev {
 public:
  template <class T>
  using array_t = std::array<T, M + 1>;

  template <class T>
  class Evaluator {
   public:
    Evaluator(const AIM::Grid &grid, const std::vector<QuantumDot> &dots)
        : eval_table{evaluation_table(grid, dots)},
          idx_table{index_table(grid, dots)},
          field_table(boost::extents[dots.size()][3])
    {
    }

    const auto &interpolate(const int time_idx,
                            const boost::multi_array<T, 6> &coef)
    {
      using mVec3_t = Eigen::Map<Eigen::Matrix<T, 3, 1>>;
      using mVec3_t_const = Eigen::Map<const Eigen::Matrix<T, 3, 1>>;
      std::fill(field_table.data(),
                field_table.data() + field_table.num_elements(), 0.0);

      enum DIMENSION { X, Y, Z };

      for(int n = 0; n < static_cast<int>(idx_table.size()); ++n) {
        mVec3_t vec(&field_table[n][0]);
        for(int i = 0; i < M + 1; ++i) {
          for(int j = 0; j < M + 1; ++j) {
            for(int k = 0; k < M + 1; ++k) {
              // clang-format off
              mVec3_t_const c(&coef[time_idx][idx_table[n]][i][j][k][0]);
              vec += Eigen::Matrix<T, 3, 1>(
                c[X] * eval_table[n][i][X][2] * eval_table[n][j][Y][0] * eval_table[n][k][Z][0] +
                c[Y] * eval_table[n][i][X][1] * eval_table[n][j][Y][1] * eval_table[n][k][Z][0] +
                c[Z] * eval_table[n][i][X][1] * eval_table[n][j][Y][0] * eval_table[n][k][Z][1],

                c[X] * eval_table[n][i][X][1] * eval_table[n][j][Y][1] * eval_table[n][k][Z][0] +
                c[Y] * eval_table[n][i][X][0] * eval_table[n][j][Y][2] * eval_table[n][k][Z][0] +
                c[Z] * eval_table[n][i][X][0] * eval_table[n][j][Y][1] * eval_table[n][k][Z][1],

                c[X] * eval_table[n][i][X][1] * eval_table[n][j][Y][0] * eval_table[n][k][Z][1] +
                c[Y] * eval_table[n][i][X][0] * eval_table[n][j][Y][1] * eval_table[n][k][Z][1] +
                c[Z] * eval_table[n][i][X][0] * eval_table[n][j][Y][0] * eval_table[n][k][Z][2]
              );
              // If you squint real hard, you'll recognize this as del del acting on the
              // tensor polnomial that fits a vector field.
              // clang-format on
            }
          }
        }
      }

      return field_table;
    }

    static boost::multi_array<double, 4> evaluation_table(
        const AIM::Grid &grid, const std::vector<QuantumDot> &dots)
    {
      using Math::Chebyshev::T;
      using Math::Chebyshev::T_d1;
      using Math::Chebyshev::T_d2;

      const Eigen::Array3d s = 2 / grid.spacing(),
                           s2 = 4 / grid.spacing().pow(2);

      std::array<int, 4> shape{{static_cast<int>(dots.size()), M + 1, 3, 3}};
      boost::multi_array<double, 4> poly_table(shape);

      for(size_t dot = 0; dot < dots.size(); ++dot) {
        Eigen::Vector3d relative_r =
            2 * (dots[dot].position().array() / grid.spacing() -
                 grid.grid_coordinate(dots[dot].position())
                     .array()
                     .cast<double>()) -
            1;
        // Effectively 2 * (r/h - floor(r/h)) - 1

        for(int n = 0; n < M + 1; ++n) {
          //          ┌────────── particle index
          //          │   ┌────── Chebyshev index (grid is outer product)
          //          │   │  ┌─── dimension (x = 0, y = 1, z = 2)
          //          │   │  │  ┌ derivative order
          //          ┴   ┴  ┴  ┴
          poly_table[dot][n][0][0] = T(n, relative_r[0]);
          poly_table[dot][n][0][1] = T_d1(n, relative_r[0]) * s[0];
          poly_table[dot][n][0][2] = T_d2(n, relative_r[0]) * s2[0];

          poly_table[dot][n][1][0] = T(n, relative_r[1]);
          poly_table[dot][n][1][1] = T_d1(n, relative_r[1]) * s[1];
          poly_table[dot][n][1][2] = T_d2(n, relative_r[1]) * s2[1];

          poly_table[dot][n][2][0] = T(n, relative_r[2]);
          poly_table[dot][n][2][1] = T_d1(n, relative_r[2]) * s[2];
          poly_table[dot][n][2][2] = T_d2(n, relative_r[2]) * s2[2];
        }
      }

      return poly_table;
    }

    static std::vector<int> index_table(const AIM::Grid &grid,
                                        const std::vector<QuantumDot> &dots)
    {
      std::vector<int> indices(dots.size());
      for(int i = 0; i < static_cast<int>(dots.size()); ++i) {
        indices[i] = grid.associated_grid_index(dots[i].position());
      }

      return indices;
    }

   private:
    boost::multi_array<double, 4> eval_table;
    std::vector<int> idx_table;
    boost::multi_array<T, 2> field_table;
  };

  Chebyshev()
      : alphas_(alphas()), roots_(roots()), poly_samples_(polynomial_samples())
  {
  }

  template <typename T>
  void fill_coefficients_tensor(const int num_boxes,
                                const T *const eval,
                                T *const coef)
  {
    constexpr int32_t size{M + 1};
    constexpr double norm{1.0 / (size * size * size)};

    const double *const chebS = &poly_samples_[0][0];

    // Generated by the Tensor Algebra Compiler (tensor-compiler.org)
    // coef(box, i, j, k, dim) =
    //    norm * alphas(i) * alphas(j) * alphas(k) *
    //    chebS(i, lambda) * chebS(j, mu) * chebS(k, nu) *
    //    eval(box, lambda, mu, nu, dim)
    for(int32_t pcoef = 0; pcoef < (num_boxes * size * size * size * 3);
        pcoef++) {
      coef[pcoef] = 0;
    }
    for(int32_t boxeval = 0; boxeval < num_boxes; boxeval++) {
      double tbox = norm;
      for(int32_t ialphas = 0; ialphas < size; ialphas++) {
        int32_t pcoef2 = boxeval * size + ialphas;
        double ti = tbox * alphas_[ialphas];
        for(int32_t lambdachebS = 0; lambdachebS < size; lambdachebS++) {
          int32_t pchebS2 = ialphas * size + lambdachebS;
          int32_t peval2 = boxeval * size + lambdachebS;
          double tlambda = ti;
          double tlambda0 = chebS[pchebS2];
          for(int32_t jalphas = 0; jalphas < size; jalphas++) {
            int32_t pcoef3 = pcoef2 * size + jalphas;
            double tj = tlambda * alphas_[jalphas];
            double tj0 = tlambda0;
            for(int32_t muchebS = 0; muchebS < size; muchebS++) {
              int32_t pchebS20 = jalphas * size + muchebS;
              int32_t peval3 = peval2 * size + muchebS;
              double tmu = tj;
              double tmu0 = tj0;
              double tmu1 = chebS[pchebS20];
              for(int32_t kalphas = 0; kalphas < size; kalphas++) {
                int32_t pcoef4 = pcoef3 * size + kalphas;
                double tk = tmu * alphas_[kalphas] * tmu0 * tmu1;
                for(int32_t nuchebS = 0; nuchebS < size; nuchebS++) {
                  int32_t pchebS21 = kalphas * size + nuchebS;
                  int32_t peval4 = peval3 * size + nuchebS;
                  double tnu = tk * chebS[pchebS21];
                  for(int32_t dimeval = 0; dimeval < 3; dimeval++) {
                    int32_t peval5 = peval4 * 3 + dimeval;
                    int32_t pcoef5 = pcoef4 * 3 + dimeval;
                    coef[pcoef5] = coef[pcoef5] + tnu * eval[peval5];
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  static array_t<int> alphas()
  {
    array_t<int> a;
    a.fill(2);
    a[0] = 1;

    return a;
  }

  static array_t<double> roots()
  {
    array_t<double> xs;
    for(int i = 0; i <= M; ++i) {
      xs[i] = -1 * std::cos(M_PI * (i + 0.5) / (M + 1));
    }

    return xs;
  }

  static array_t<array_t<double>> polynomial_samples()
  {
    array_t<array_t<double>> table;

    const auto xk = roots();

    for(int k = 0; k <= M; ++k) {
      table[0][k] = 1;
      if(M > 0) table[1][k] = xk[k];
    }

    for(int p = 2; p <= M; ++p) {
      for(int k = 0; k <= M; ++k) {
        table[p][k] = 2 * xk[k] * table[p - 1][k] - table[p - 2][k];
      }
    }

    return table;
  }

 private:
  array_t<int> alphas_;
  array_t<double> roots_;
  array_t<array_t<double>> poly_samples_;
};

#endif
