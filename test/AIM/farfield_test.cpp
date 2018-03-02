#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iomanip>

#include "interactions/AIM/farfield.h"

BOOST_AUTO_TEST_SUITE(AIM)

BOOST_AUTO_TEST_SUITE(FARFIELD)

BOOST_AUTO_TEST_CASE(CONSTRUCTION)
{
  auto dots = std::make_shared<DotVector>();
  dots->push_back(QuantumDot({0.1, 0.1, 0.1}));
  dots->push_back(QuantumDot({0.1, 0.1, 10.1}));

  auto history =
      std::make_shared<Integrator::History<Eigen::Vector2cd>>(2, 10, 1024);

  AIM::Grid grid({1, 1, 1}, 1, *dots);

  AIM::Expansions::LeastSquaresExpansionSolver LSE(grid);
  auto expansion_table = LSE.table(*dots);
  auto cheb_expansion_table = LSE.chebyshev_lambda_weights(
      Math::Chebyshev::normalized_points(AIM::chebyshev_order));

  AIM::Farfield(dots, history, 4, 1, 1, grid, expansion_table,
                AIM::normalization::unit, cheb_expansion_table);
}

BOOST_AUTO_TEST_SUITE_END()  // FARFIELD

BOOST_AUTO_TEST_SUITE_END()  // AIM
