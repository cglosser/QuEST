#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iomanip>

#include "interactions/AIM/farfield.h"
#include "interactions/AIM/nearfield.h"

BOOST_AUTO_TEST_SUITE(AIM)

struct PARAMETERS {
  using Hist_t = Integrator::History<Eigen::Vector2cd>;
  using LSE = AIM::Expansions::LeastSquaresExpansionSolver;

  std::vector<Eigen::Vector3d> default_pos;

  int n_steps, cheb_order;

  std::shared_ptr<DotVector> dots;
  std::shared_ptr<Hist_t> history;

  PARAMETERS()
      : default_pos{{0.57, 0.28, 0.04}, {0.76, 0.48, 0.41}, {0.78, 0.56, 0.24},
                    {0.42, 0.66, 0.70}, {0.07, 0.77, 0.77}, {0.68, 0.99, 0.45}},
        n_steps{1024},
        cheb_order{3},
        dots{std::make_shared<DotVector>()},
        history{std::make_shared<Hist_t>(default_pos.size(), 10, n_steps, 1)}
  {
    history->fill(Eigen::Vector2cd::Zero());
    for(int t = -10; t < n_steps; ++t) {
      history->array_[0][t][0] = Eigen::Vector2cd(
          0, Math::gaussian((t - n_steps / 2.0) / (n_steps / 12.0)));
    }
  }

  AIM::Grid make_grid() { return AIM::Grid({1, 1, 1}, 1, *dots); }
  auto make_interactions()
  {
    auto grid{make_grid()};
    LSE lse(grid);
    auto expansion_table = lse.table(*dots);
    auto cheb_table = lse.chebyshev_lambda_weights(
        Math::Chebyshev::normalized_points(cheb_order));

    return std::make_pair(
        std::make_unique<AIM::Nearfield>(dots, history, 4, 100, 1, 1, grid,
                                         expansion_table,
                                         AIM::Normalization::unit, cheb_table),
        std::make_unique<AIM::Farfield>(dots, history, 4, 1, 1, grid,
                                        expansion_table,
                                        AIM::Normalization::unit, cheb_table));
  }
};

BOOST_FIXTURE_TEST_SUITE(FIELDS, PARAMETERS)

BOOST_AUTO_TEST_CASE(SAME_BOX)
{
  for(const auto &r : default_pos) {
    dots->emplace_back(r, Eigen::Vector3d(0, 0, 1));
  }

  auto inter{make_interactions()};

  for(int t = 0; t < n_steps; ++t) {
    auto eval = inter.first->evaluate(t) - inter.second->evaluate(t);
    BOOST_CHECK_SMALL(eval.matrix().norm(), 1e-12);
  }
}

BOOST_AUTO_TEST_SUITE_END()  // FIELDS

BOOST_AUTO_TEST_SUITE_END()  // AIM
