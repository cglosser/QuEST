#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iomanip>

#include "interactions/AIM/farfield.h"

BOOST_AUTO_TEST_SUITE(AIM)

BOOST_AUTO_TEST_SUITE(FARFIELD)

BOOST_AUTO_TEST_CASE(CONSTRUCTION)
{
  auto dots = std::make_shared<DotVector>();
  dots->push_back(QuantumDot({0.1, 0.1, 0.1}, {0, 0, 1}));
  //dots->push_back(QuantumDot({0.1, 0.1, 3.1}, {0, 0, 1}));

  auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(
      dots->size(), 10, 1024);
  for(int t = -10; t < 1024; ++t) {
    history->array_[0][t][0] =
        Eigen::Vector2cd(0, Math::gaussian((t - 512) / (1024.0 / 12)));
    //history->array_[1][t][0] =
        //Eigen::Vector2cd(0, -0.5 * Math::gaussian((t - 512) / (1024.0 / 12)));
  }

  AIM::Grid grid({1, 1, 1}, 1, *dots);

  AIM::Expansions::LeastSquaresExpansionSolver LSE(grid);
  auto expansion_table = LSE.table(*dots);
  auto cheb_expansion_table = LSE.chebyshev_lambda_weights(
      Math::Chebyshev::normalized_points(AIM::chebyshev_order));

  AIM::Farfield ff(dots, history, 4, 1, 1, grid, expansion_table,
                   AIM::normalization::unit, cheb_expansion_table);

  std::cout << std::setprecision(14);
  for(int t = 0; t < 1024; ++t) {
    std::cout << ff.evaluate(t).transpose() << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()  // FARFIELD

BOOST_AUTO_TEST_SUITE_END()  // AIM
