#ifndef PROJECTORS_H
#define PROJECTORS_H

namespace Projector {
  enum DIMENSION { X, Y, Z };

  template <typename T>
  auto Potential(const int t,
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

  template <typename T>
  auto GradDiv(const int t,
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

  template <typename T>
  class TimeDerivative {
   public:
    TimeDerivative(const int history_length, const double dt = 1)
        : history_length(history_length),
          dt_coefs({{25.0 / 12, -4.0, 3.0, -4.0 / 3, 1.0 / 4}})
    {
      for(auto &c : dt_coefs) c /= dt;
    }

    auto operator()(const int t,
                    const int n,
                    const int i,
                    const int j,
                    const int k,
                    const boost::multi_array<T, 6> &coef,
                    const boost::multi_array<double, 4> &eval) const
    {
      Eigen::Matrix<T, 3, 1> total_field;
      total_field.setZero();

      double Ts = eval[n][i][X][0] * eval[n][j][Y][0] * eval[n][k][Z][0];

      for(int h = 0; h < static_cast<int>(dt_coefs.size()); ++h) {
        int w = std::max(t - h, 0) % history_length;
        Eigen::Map<const Eigen::Matrix<T, 3, 1>> c(&coef[w][n][i][j][k][0]);
        total_field += dt_coefs[h] * c * Ts;
      }
      return total_field;
    }

   private:
    int history_length;
    std::array<double, 5> dt_coefs;
  };
}

#endif
