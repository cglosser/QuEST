#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iomanip>

#include "interactions/AIM/aim_interaction.h"

BOOST_AUTO_TEST_SUITE(AIM)

BOOST_AUTO_TEST_SUITE(FARFIELD)

BOOST_AUTO_TEST_CASE(CONSTRUCTOR)
{
  auto dots = std::make_shared<DotVector>();
  dots->push_back(QuantumDot({0, 0, 0}, {0, 0, 1}));

  using hist_t = Integrator::History<Eigen::Vector2cd>;
  auto hist = std::make_shared<hist_t>(1, 10, 1024, 1);

  auto grid = std::make_shared<AIM::Grid>(Eigen::Array3d(1, 1, 1), 1, *dots);

  using LSE = AIM::Expansions::LeastSquaresExpansionSolver;
  auto expansion_table = std::make_shared<AIM::Expansions::ExpansionTable>(
      LSE::get_expansions(1, *grid, *dots));

  AIM::Expansions::Retardation potential(grid->max_transit_steps(1, 1) + 4);

  AIM::Farfield ff(dots, hist, 4, 1, 1, grid, expansion_table, potential,
                   AIM::Normalization::unit);
}

BOOST_AUTO_TEST_SUITE_END()  // FARFIELD

BOOST_AUTO_TEST_SUITE_END()  // AIM
