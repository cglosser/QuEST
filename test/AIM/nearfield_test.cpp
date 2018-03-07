#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iomanip>

#include "interactions/AIM/nearfield.h"

BOOST_AUTO_TEST_SUITE(AIM)

BOOST_AUTO_TEST_SUITE(NEARFIELD)

BOOST_AUTO_TEST_CASE(CONSTRUCTION)
{
  auto dots = std::make_shared<DotVector>();
  dots->push_back(QuantumDot(
      {0.3086582838174551, 0.3086582838174551, 0.3086582838174551}, {0, 0, 1}));
  dots->push_back(
      QuantumDot({0.961939766255643, 0.961939766255643, 10 + 0.961939766255643},
                 {0, 0, 1}));

  auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(
      dots->size(), 10, 1024);
  for(int t = -10; t < 1024; ++t) {
    history->array_[0][t][0] =
        Eigen::Vector2cd(0, Math::gaussian((t - 512) / (1024.0 / 12)));
  }

  AIM::Grid grid({1, 1, 1}, 1, *dots);

  AIM::Expansions::LeastSquaresExpansionSolver LSE(grid);
  auto expansion_table = LSE.table(*dots);
  auto cheb_expansion_table = LSE.chebyshev_lambda_weights(
      Math::Chebyshev::normalized_points(AIM::chebyshev_order));
  Projector::Potential<cmplx> proj(grid.max_transit_steps(1, 1) + 4);

  AIM::Nearfield nf(dots, history, 4, 100, 1, 1, grid, expansion_table,
                    AIM::Normalization::unit, cheb_expansion_table, proj);

  std::cout.width(14);
  std::cout << std::setprecision(14);
  for(int t = 0; t < 1024; ++t) {
    std::cout << nf.evaluate(t).transpose() << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()  // FARFIELD

BOOST_AUTO_TEST_SUITE_END()  // AIM
