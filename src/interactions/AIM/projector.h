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
  class Derivative : public ProjectorBase<num_t> {
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
      Eigen::Map<const Eigen::Array<num_t, 3, 1>> c(
          &coef[this->wrap_step(t)][box][i][j][k][0]);

      double deriv = eval[n][i][X][0] * eval[n][j][Y][0] * eval[n][k][Z][1];
      return deriv * c;
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

  template <typename num_t>
  class EFIE : public ProjectorBase<num_t> {
   public:
    EFIE(const int history_length, const double c = 1, const double dt = 1)
        : ProjectorBase<num_t>(history_length),
          c_{c},
          dt_sq_coefs{{35.0 / 12, -26.0 / 3, 19.0 / 2, -14.0 / 3, 11.0 / 12}}
    {
      for(auto &x : dt_sq_coefs) x /= std::pow(dt, 2);
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
      double Ts = eval[n][i][X][0] * eval[n][j][Y][0] * eval[n][k][Z][0];
      Eigen::Matrix<num_t, 3, 1> dt_sq = Eigen::Matrix<num_t, 3, 1>::Zero();

      for(int h = 0; h < static_cast<int>(dt_sq_coefs.size()); ++h) {
        int w = this->wrap_step(std::max(t - h, 0));
        Eigen::Map<const Eigen::Matrix<num_t, 3, 1>> c(
            &coef[w][box][i][j][k][0]);
        dt_sq += dt_sq_coefs[h] * c * Ts;
      }

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

      return dt_sq -
             std::pow(c_, 2) * m *
                 Eigen::Map<const Eigen::Matrix<num_t, 3, 1>>(
                     &coef[this->wrap_step(t)][box][i][j][k][0]);
    }

   private:
    double c_;
    std::array<double, 5> dt_sq_coefs;
  };

  class RotatingEFIE : public ProjectorBase<cmplx> {
   public:
    RotatingEFIE(const int history_length,
                 const double c = 1,
                 const double omega = 0,
                 const double dt = 1)
        : ProjectorBase<cmplx>(history_length),
          c_{c},
          omega_{omega},
          dt_coefs{{25.0 / 12, -4.0, 3, -4.0 / 3, 1.0 / 4}},
          dt_sq_coefs{{35.0 / 12, -26.0 / 3, 19.0 / 2, -14.0 / 3, 11.0 / 12}}
    {
      for(auto &x : dt_coefs) x /= dt;
      for(auto &x : dt_sq_coefs) x /= std::pow(dt, 2);
    }

    Eigen::Matrix<cmplx, 3, 1> operator()(
        const int t,
        const int n,
        const int box,
        const int i,
        const int j,
        const int k,
        const boost::multi_array<cmplx, 6> &coef,
        const boost::multi_array<double, 4> &eval) const
    {
      Eigen::Matrix3cd deriv_t;
      deriv_t.setZero();

      Eigen::Map<const Eigen::Vector3cd> current(
          &coef[this->wrap_step(t)][box][i][j][k][0]);

      double Ts = eval[n][i][X][0] * eval[n][j][Y][0] * eval[n][k][Z][0];
      deriv_t.col(0) = Ts * current;

      for(int h = 0; h < static_cast<int>(dt_coefs.size()); ++h) {
        int w = this->wrap_step(std::max(t - h, 0));
        Eigen::Map<const Eigen::Vector3cd> c(&coef[w][box][i][j][k][0]);
        deriv_t.col(1) += dt_coefs[h] * c * Ts;
        deriv_t.col(2) += dt_sq_coefs[h] * c * Ts;
      }

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

      return (deriv_t.col(2) + 2.0 * iu * omega_ * deriv_t.col(1) -
              std::pow(omega_, 2) * deriv_t.col(0)) -
             std::pow(c_, 2) * m * current;
    }

   private:
    double c_, omega_;
    std::array<double, 5> dt_coefs, dt_sq_coefs;
  };
}

#endif
