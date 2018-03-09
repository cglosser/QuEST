#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iomanip>

#include "interactions/AIM/farfield.h"
#include "interactions/AIM/nearfield.h"

BOOST_AUTO_TEST_SUITE(AIM)

struct PARAMETERS {
  using Hist_t = Integrator::History<Eigen::Vector2cd>;
  using LSE = AIM::Expansions::LeastSquaresExpansionSolver;

  int n_pts, n_steps, cheb_order;
  double c, dt;

  std::shared_ptr<DotVector> dots;
  std::shared_ptr<Hist_t> history;

  PARAMETERS(int n_pts, int n_steps)
      : n_pts{n_pts},
        n_steps{n_steps},
        cheb_order{AIM::chebyshev_order},
        c{1},
        dt{1},
        dots{std::make_shared<DotVector>()},
        history{std::make_shared<Hist_t>(n_pts, 10, n_steps, 1)}
  {
    std::cout << std::setprecision(17) << std::scientific;
  }
};

BOOST_AUTO_TEST_SUITE(ANALYTIC_COMPARISON)

class Gaussian {
 public:
  Gaussian(double mu, double sigma) : mu_{mu}, sigma_{sigma} {};

  double operator()(int n, double x)
  {
    switch(n) {
      case 0: return dt0(x);
      case 1: return dt1(x);
      case 2: return dt2(x);
      default: return 0;
    }
  }

 private:
  double dt0(double x) const { return Math::gaussian((x - mu_) / sigma_); }
  double dt1(double x) const
  {
    double arg = (x - mu_) / sigma_;
    return -1 * arg * Math::gaussian(arg) / sigma_;
  }
  double dt2(double x) const
  {
    double arg = (x - mu_) / sigma_;
    return (std::pow(arg, 2) - 1) * Math::gaussian(arg) / std::pow(sigma_, 2);
  }

  double mu_, sigma_;
};

struct ANALYTIC_PARAMETERS : public PARAMETERS {
  Gaussian gaussian;
  Eigen::Vector3d dr;
  double dr_norm;

  ANALYTIC_PARAMETERS()
      : PARAMETERS(2, 1024), gaussian(n_steps / 2.0, n_steps / 12.0)
  {
    Eigen::Vector3d z_hat(0, 0, 1);
    dots->emplace_back(Eigen::Vector3d(0.1, 0.1, 0.1), z_hat);
    dots->emplace_back(Eigen::Vector3d(0.9, 0.9, 5.9), z_hat);
    dr = dots->at(1).position() - dots->at(0).position();
    dr_norm = dr.norm();

    history->fill(Eigen::Vector2cd::Zero());
    for(int t = -10; t < n_steps; ++t) {
      history->array_[0][t][0] = Eigen::Vector2cd(0, gaussian(0, t * dt));
    }
  }

  AIM::Grid make_grid() { return AIM::Grid({0.7, 0.7, 0.7}, 2, *dots); }
  template <typename P>
  auto make_interaction(AIM::Normalization::SpatialNorm norm)
  {
    auto grid{make_grid()};
    LSE lse(grid);
    auto expansion_table =
        std::make_shared<AIM::Expansions::ExpansionTable>(lse.table(*dots));
    auto cheb_table = lse.chebyshev_lambda_weights(
        Math::Chebyshev::normalized_points(cheb_order));
    int interp_order = 4;

    P proj(grid.max_transit_steps(c, dt) + interp_order);

    return std::make_unique<AIM::Nearfield>(dots, history, interp_order, 100, c,
                                            dt, grid, expansion_table, norm,
                                            cheb_table, proj);
  }

  template <typename T>
  double test_analytic(T &interaction, std::function<cmplx(double)> solution)
  {
    double max_err = 0;
    for(int t = 0; t < n_steps; ++t) {
      auto x = interaction->evaluate(t);
      double err = std::abs(solution(t) - x(1));
      if(t > n_steps / 10) max_err = std::max(max_err, err);
    }

    BOOST_TEST_MESSAGE("Maximum error: " << max_err);
    return max_err;
  }
};

BOOST_FIXTURE_TEST_CASE(RETARDATION, ANALYTIC_PARAMETERS)
{
  auto solution = [&](double t) -> cmplx {
    return cmplx(gaussian(0, t - dr_norm / c), 0);
  };

  auto nf =
      make_interaction<Projector::Potential<cmplx>>(AIM::Normalization::unit);

  double err = test_analytic(nf, solution);
  BOOST_CHECK_SMALL(err, 1.0);
}

BOOST_FIXTURE_TEST_CASE(TIME_DERIVATIVE, ANALYTIC_PARAMETERS)
{
  auto solution = [&](double t) -> cmplx {
    return cmplx(gaussian(1, t - dr_norm / c), 0);
  };

  auto nf = make_interaction<Projector::TimeDerivative<cmplx>>(
      AIM::Normalization::unit);

  double err = test_analytic(nf, solution);
  BOOST_CHECK_SMALL(err, 1.0);
}

BOOST_FIXTURE_TEST_CASE(GRAD_DIV, ANALYTIC_PARAMETERS)
{
  auto solution = [&](double t) -> cmplx {
    double arg = t - dr_norm / c;
    double val =
        (c * (std::pow(dr(2), 2) - std::pow(dr_norm, 2)) * gaussian(1, arg) +
         dr_norm * std::pow(dr(2), 2) * gaussian(2, arg)) /
        (std::pow(c, 2) * std::pow(dr_norm, 3));
    return {val, 0};
  };

  auto nf =
      make_interaction<Projector::GradDiv<cmplx>>(AIM::Normalization::unit);

  double err = test_analytic(nf, solution);
  BOOST_CHECK_SMALL(err, 1.0);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(EQUIVALENT_NEAR_AND_FAR_FIELDS)

struct EQUIVALENCE_BASE : public PARAMETERS {
  struct FieldPair {
    std::unique_ptr<AIM::Nearfield> nf;
    std::unique_ptr<AIM::Farfield> ff;
  };

  std::vector<Eigen::Vector3d> default_pos;

  EQUIVALENCE_BASE()
      : PARAMETERS(6, 256),
        default_pos{{0.57, 0.28, 0.04}, {0.76, 0.48, 0.41}, {0.78, 0.56, 0.24},
                    {0.42, 0.66, 0.70}, {0.07, 0.77, 0.77}, {0.68, 0.99, 0.45}}
  {
    history->fill(Eigen::Vector2cd::Zero());
    for(int i = 0; i < n_pts; ++i) {
      for(int t = -10; t < n_steps; ++t) {
        history->array_[i][t][0] = Eigen::Vector2cd(
            0, Math::gaussian((t - n_steps / 2.0) / (n_steps / 12.0)));
      }
    }
  }

  AIM::Grid make_grid() { return AIM::Grid({1, 1, 1}, 1, *dots); }
  void test_equivalence(const FieldPair &p)
  {
    for(int t = 0; t < n_steps; ++t) {
      auto eval = p.nf->evaluate(t) - p.ff->evaluate(t);
      BOOST_CHECK_SMALL(eval.matrix().norm(), 1e-12);
    }
  }
};

struct RETARDATION_PARAMETERS : public EQUIVALENCE_BASE {
  auto make_interactions()
  {
    auto grid{make_grid()};
    LSE lse(grid);
    auto expansion_table =
        std::make_shared<AIM::Expansions::ExpansionTable>(lse.table(*dots));
    auto cheb_table = lse.chebyshev_lambda_weights(
        Math::Chebyshev::normalized_points(cheb_order));
    int interp_order = 4;

    Projector::Potential<cmplx> potential(grid.max_transit_steps(c, dt) +
                                          interp_order);

    return FieldPair{
        std::make_unique<AIM::Nearfield>(
            dots, history, interp_order, 100, c, dt, grid, expansion_table,
            AIM::Normalization::unit, cheb_table, potential),
        std::make_unique<AIM::Farfield>(
            dots, history, interp_order, c, dt, grid, expansion_table,
            AIM::Normalization::unit, cheb_table, potential)};
  }
};

BOOST_FIXTURE_TEST_SUITE(RETARDATION, RETARDATION_PARAMETERS)

BOOST_AUTO_TEST_CASE(SAME_BOX)
{
  for(const auto &r : default_pos) {
    dots->emplace_back(r, Eigen::Vector3d(0, 0, 1));
  }

  FieldPair fieldPair{make_interactions()};
  test_equivalence(fieldPair);
}

BOOST_AUTO_TEST_CASE(ADJACENT_BOX)
{
  Eigen::Vector3d z_hat(0, 0, 1);
  for(int i = 0; i < n_pts; ++i) {
    dots->emplace_back(default_pos[i] + (i % 2) * z_hat, z_hat);
  }

  FieldPair fieldPair{make_interactions()};
  test_equivalence(fieldPair);
}

BOOST_AUTO_TEST_CASE(DISTANT_BOX)
{
  Eigen::Vector3d z_hat(0, 0, 1);
  for(int i = 0; i < n_pts; ++i) {
    dots->emplace_back(default_pos[i] + (i % 2) * 10 * z_hat, z_hat);
  }

  FieldPair fieldPair{make_interactions()};
  test_equivalence(fieldPair);
}

BOOST_AUTO_TEST_SUITE_END()  // RETARDATION

struct LAPLACE_PARAMETERS : public EQUIVALENCE_BASE {
  auto make_interactions()
  {
    auto grid{make_grid()};
    LSE lse(grid);
    auto expansion_table =
        std::make_shared<AIM::Expansions::ExpansionTable>(lse.table(*dots));
    auto cheb_table = lse.chebyshev_lambda_weights(
        Math::Chebyshev::normalized_points(cheb_order));
    int interp_order = 4;

    Projector::Potential<cmplx> potential(grid.max_transit_steps(c, dt) +
                                          interp_order);

    return FieldPair{std::make_unique<AIM::Nearfield>(
                         dots, history, 4, 100, 1, 1, grid, expansion_table,
                         AIM::Normalization::Laplace(), cheb_table, potential),
                     std::make_unique<AIM::Farfield>(
                         dots, history, 4, 1, 1, grid, expansion_table,
                         AIM::Normalization::Laplace(), cheb_table, potential)};
  }
};

BOOST_FIXTURE_TEST_SUITE(LAPLACE, LAPLACE_PARAMETERS)

BOOST_AUTO_TEST_CASE(SAME_BOX)
{
  for(const auto &r : default_pos) {
    dots->emplace_back(r, Eigen::Vector3d(0, 0, 1));
  }

  FieldPair fieldPair{make_interactions()};
  test_equivalence(fieldPair);
}

BOOST_AUTO_TEST_CASE(ADJACENT_BOX)
{
  Eigen::Vector3d z_hat(0, 0, 1);
  for(int i = 0; i < n_pts; ++i) {
    dots->emplace_back(default_pos[i] + (i % 2) * z_hat, z_hat);
  }

  FieldPair fieldPair{make_interactions()};
  test_equivalence(fieldPair);
}

BOOST_AUTO_TEST_CASE(DISTANT_BOX)
{
  Eigen::Vector3d z_hat(0, 0, 1);
  for(int i = 0; i < n_pts; ++i) {
    dots->emplace_back(default_pos[i] + (i % 2) * 10 * z_hat, z_hat);
  }

  FieldPair fieldPair{make_interactions()};
  test_equivalence(fieldPair);
}

BOOST_AUTO_TEST_SUITE_END()  // LAPLACE

BOOST_AUTO_TEST_SUITE_END()  // EQUIVALENT_NEAR_AND_FAR_FIELDS

BOOST_AUTO_TEST_SUITE_END()  // AIM
