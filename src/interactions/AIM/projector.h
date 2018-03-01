#ifndef PROJECTORS_H
#define PROJECTORS_H

namespace Projector {
  enum DIMENSION { X, Y, Z };

  template <class T>
  auto potential(const int t,
                 const int n,
                 const int i,
                 const int j,
                 const int k,
                 const boost::multi_array<T, 6> &coef,
                 const boost::multi_array<double, 4> &eval)
  {
    Eigen::Map<const Eigen::Matrix<T, 3, 1>> c(&coef[t][n][i][j][k][0]);
    double Ts = eval[n][i][X][0] * eval[n][j][Y][0] * eval[n][k][Z][0];
    return c * Ts;
  }

  template <class T>
  auto grad_div(const int t,
                const int n,
                const int i,
                const int j,
                const int k,
                const boost::multi_array<T, 6> &coef,
                const boost::multi_array<double, 4> &eval)
  {
    Eigen::Map<const Eigen::Matrix<T, 3, 1>> c(&coef[t][n][i][j][k][0]);
    Eigen::Matrix3d m;
    m << eval[n][i][X][2] * eval[n][j][Y][0] * eval[n][k][Z][0],
        eval[n][i][X][1] * eval[n][j][Y][1] * eval[n][k][Z][0],
        eval[n][i][X][1] * eval[n][j][Y][1] * eval[n][k][Z][0],

        eval[n][i][X][1] * eval[n][j][Y][1] * eval[n][k][Z][0],
        eval[n][i][X][0] * eval[n][j][Y][2] * eval[n][k][Z][0],
        eval[n][i][X][0] * eval[n][j][Y][1] * eval[n][k][Z][1],

        eval[n][i][X][1] * eval[n][j][Y][0] * eval[n][k][Z][1],
        eval[n][i][X][0] * eval[n][j][Y][1] * eval[n][k][Z][1],
        eval[n][i][X][0] * eval[n][j][Y][0] * eval[n][k][Z][2];

    return (m * c).eval();
  }
}

#endif
