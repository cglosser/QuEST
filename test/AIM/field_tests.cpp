#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iomanip>

#include "interactions/AIM/farfield.h"
#include "interactions/AIM/nearfield.h"

BOOST_AUTO_TEST_SUITE(AIM)

struct PARAMETERS {
  using Hist_t = Integrator::History<Eigen::Vector2cd>;
  using LSE = AIM::Expansions::LeastSquaresExpansionSolver;

  int n_pts, n_steps;
  double c, dt;

  std::shared_ptr<DotVector> dots;
  std::shared_ptr<Hist_t> history;

  PARAMETERS(int n_pts, int n_steps)
      : n_pts{n_pts},
        n_steps{n_steps},
        c{1},
        dt{1},
        dots{std::make_shared<DotVector>()},
        history{std::make_shared<Hist_t>(n_pts, 10, n_steps, 1)}
  {
    std::cout << std::setprecision(17) << std::scientific;
  }
};

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
            0, gaussian((t - n_steps / 2.0) / (n_steps / 12.0)));
      }
    }
  }

  auto make_grid() { return AIM::Grid({1, 1, 1}, 1, *dots); }
  void test_equivalence(const FieldPair &p)
  {
    for(int t = 0; t < n_steps; ++t) {
      auto eval = p.nf->evaluate(t) - p.ff->evaluate(t);
      if(t > n_steps / 10) BOOST_CHECK_SMALL(eval.matrix().norm(), 1e-12);
    }
  }
};

struct RETARDATION_PARAMETERS : public EQUIVALENCE_BASE {
  auto make_interactions()
  {
    auto grid = std::make_shared<AIM::Grid>(make_grid());
    auto expansion_table = std::make_shared<AIM::Expansions::ExpansionTable>(
        LSE::get_expansions(1, *grid, *dots));
    auto nf_pairs = std::make_shared<std::vector<AIM::Grid::ipair_t>>(
        grid->nearfield_point_pairs(100, *dots));

    int interp_order = 4;

    AIM::Expansions::Retardation potential(grid->max_transit_steps(c, dt) +
                                           interp_order);

    return FieldPair{
        std::make_unique<AIM::Nearfield>(dots, history, interp_order, c, dt,
                                         grid, expansion_table, potential,
                                         AIM::Normalization::unit, nf_pairs),
        std::make_unique<AIM::Farfield>(dots, history, interp_order, c, dt,
                                        grid, expansion_table, potential,
                                        AIM::Normalization::unit)};
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
    auto grid{std::make_shared<AIM::Grid>(make_grid())};
    auto expansion_table = std::make_shared<AIM::Expansions::ExpansionTable>(
        LSE::get_expansions(1, *grid, *dots));
    auto nf_pairs = std::make_shared<std::vector<AIM::Grid::ipair_t>>(
        grid->nearfield_point_pairs(100, *dots));
    int interp_order = 4;

    AIM::Expansions::Retardation potential(grid->max_transit_steps(c, dt) +
                                           interp_order);

    return FieldPair{std::make_unique<AIM::Nearfield>(
                         dots, history, 4, 1, 1, grid, expansion_table,
                         potential, AIM::Normalization::Laplace(), nf_pairs),
                     std::make_unique<AIM::Farfield>(
                         dots, history, 4, 1, 1, grid, expansion_table,
                         potential, AIM::Normalization::Laplace())};
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

// BOOST_AUTO_TEST_CASE(ADJACENT_BOX)
//{
// Eigen::Vector3d z_hat(0, 0, 1);
// for(int i = 0; i < n_pts; ++i) {
// dots->emplace_back(default_pos[i] + (i % 2) * z_hat, z_hat);
//}

// FieldPair fieldPair{make_interactions()};
// test_equivalence(fieldPair);
//}

// BOOST_AUTO_TEST_CASE(DISTANT_BOX)
//{
// Eigen::Vector3d z_hat(0, 0, 1);
// for(int i = 0; i < n_pts; ++i) {
// dots->emplace_back(default_pos[i] + (i % 2) * 10 * z_hat, z_hat);
//}

// FieldPair fieldPair{make_interactions()};
// test_equivalence(fieldPair);
//}

BOOST_AUTO_TEST_SUITE_END()  // LAPLACE

BOOST_AUTO_TEST_SUITE_END()  // EQUIVALENT_NEAR_AND_FAR_FIELDS

BOOST_AUTO_TEST_SUITE_END()  // AIM
