#ifndef PROJECTORS_H
#define PROJECTORS_H

namespace Projector {
  enum DIMENSION { X, Y, Z };

  template <typename num_t>
  using Projector_fn = std::function<Eigen::Matrix<num_t, 3, 1>(
      const int,
      const int,
      const int,
      const int,
      const int,
      const int,
      const boost::multi_array<num_t, 6> &,
      const boost::multi_array<double, 4> &)>;

  template <typename num_t>
  class ProjectorBase {
   public:
    ProjectorBase(int max_timestep) : max_timestep_{max_timestep} {}
   protected:
    int wrap_step(int t) const { return t % max_timestep_; }
    int max_timestep_;
  };

  template <typename num_t>
  class Potential : public ProjectorBase<num_t> {
   public:
    using ProjectorBase<num_t>::ProjectorBase;
    Eigen::Matrix<num_t, 3, 1> operator()(
        const int t,
        const int n,
        const int box,
        const int i,
        const int j,
        const int k,
        const boost::multi_array<num_t, 6> &coef,
        const boost::multi_array<double, 4> &eval)
    {
      Eigen::Map<const Eigen::Matrix<num_t, 3, 1>> c(
          &coef[this->wrap_step(t)][box][i][j][k][0]);
      double Ts = eval[n][i][X][0] * eval[n][j][Y][0] * eval[n][k][Z][0];
      return c * Ts;
    }
  };

  template <typename num_t>
  class GradDiv : public ProjectorBase<num_t> {
   public:
    using ProjectorBase<num_t>::ProjectorBase;
    Eigen::Matrix<num_t, 3, 1> operator()(
        const int t,
        const int n,
        const int box,
        const int i,
        const int j,
        const int k,
        const boost::multi_array<num_t, 6> &coef,
        const boost::multi_array<double, 4> &eval)
    {
      Eigen::Map<const Eigen::Matrix<num_t, 3, 1>> c(
          &coef[this->wrap_step(t)][box][i][j][k][0]);
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
  };

  template <typename num_t>
  class TimeDerivative : public ProjectorBase<num_t> {
   public:
    TimeDerivative(const int history_length, const double dt = 1)
        : ProjectorBase<num_t>(history_length),
          dt_coefs({{25.0 / 12, -4.0, 3.0, -4.0 / 3, 1.0 / 4}})
    {
      for(auto &c : dt_coefs) c /= dt;
    }

    Eigen::Matrix<num_t, 3, 1> operator()(
        const int t,
        const int n,
        const int box,
        const int i,
        const int j,
        const int k,
        const boost::multi_array<num_t, 6> &coef,
        const boost::multi_array<double, 4> &eval) const
    {
      Eigen::Matrix<num_t, 3, 1> total_field;
      total_field.setZero();

      double Ts = eval[n][i][X][0] * eval[n][j][Y][0] * eval[n][k][Z][0];

      for(int h = 0; h < static_cast<int>(dt_coefs.size()); ++h) {
        int w = this->wrap_step(std::max(t - h, 0));
        Eigen::Map<const Eigen::Matrix<num_t, 3, 1>> c(
            &coef[w][box][i][j][k][0]);
        total_field += dt_coefs[h] * c * Ts;
      }
      return total_field;
    }

   private:
    std::array<double, 5> dt_coefs;
  };
}

#endif
