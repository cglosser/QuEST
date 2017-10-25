#include <boost/test/unit_test.hpp>

#include "interactions/AIM/expansion.h"
#include "interactions/AIM/grid.h"

BOOST_AUTO_TEST_SUITE(EXPANSIONS)

struct ORDERS {
  int interpolation_order, expansion_order;
  ORDERS() : interpolation_order(4), expansion_order(1){};
};

BOOST_FIXTURE_TEST_CASE(ONE_POINT_EXPANSION, ORDERS)
{
  std::shared_ptr<DotVector> dots = std::make_shared<DotVector>(
      DotVector{QuantumDot(Eigen::Vector3d(0.5, 0.5, 0.5))});
  AIM::Grid grid(Eigen::Vector3d(1, 1, 1), dots, expansion_order);
  auto expansions = AIM::LeastSquaresExpansionSolver::get_expansions(
      expansion_order, grid, *dots);

  for(int i = 0; i < 8; ++i) {
    BOOST_CHECK_EQUAL(expansions[0][i].weight, 1.0 / 8);
  }
}

BOOST_FIXTURE_TEST_CASE(TWO_POINT_EXPANSION, ORDERS)
{
  std::shared_ptr<DotVector> dots = std::make_shared<DotVector>(
      DotVector{QuantumDot(Eigen::Vector3d(0.5, 0.5, 0.5)),
                QuantumDot(Eigen::Vector3d(0.5, 0.5, 10.5))});
  AIM::Grid grid(Eigen::Vector3d(1, 1, 1), dots, expansion_order);
  auto expansions = AIM::LeastSquaresExpansionSolver::get_expansions(
      expansion_order, grid, *dots);

  for(int dot = 0; dot < 2; ++dot) {
    for(int pt = 0; pt < 8; ++pt) {
      BOOST_CHECK_EQUAL(expansions[dot][pt].weight, 1.0 / 8);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()  // EXPANSIONS
