#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iomanip>

#include "interactions/AIM/aim_interaction.h"

BOOST_AUTO_TEST_SUITE(AIM)

BOOST_AUTO_TEST_SUITE(NEARFIELD)

struct PARAMETERS {
  int interpolation_order, num_steps, num_dots;
  double c, dt, total_time;

  Eigen::Array3d spacing;

  int expansion_order;  // Different orders for different test geometries

  std::shared_ptr<DotVector> dots;
  std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history;

  PARAMETERS(const int expansion_order)
      : interpolation_order(3),
        num_steps(1024),
        num_dots(2),

        c(1),
        dt(1),
        total_time(dt * num_steps),

        spacing(Eigen::Array3d(1, 1, 1) * c * dt),

        expansion_order(expansion_order),
        dots(std::make_shared<DotVector>()),
        history(std::make_shared<Integrator::History<Eigen::Vector2cd>>(
            num_dots, 10, num_steps)){};

  double src(const double t) const
  {
    double arg = (t - total_time / 2.0) / (total_time / 6.0);
    return gaussian(arg);
  }
};

struct ON_GRID_PARAMETERS : public PARAMETERS {
  ON_GRID_PARAMETERS() : PARAMETERS(1)
  {
    history->fill(Eigen::Vector2cd::Zero());
    for(int i = -10; i < num_steps; ++i) {
      history->array[0][i][0](RHO_01) = src(i * dt);
      history->array[1][i][0](RHO_01) = 0;
    }
  }
};

BOOST_FIXTURE_TEST_CASE(TWO_PC, ON_GRID_PARAMETERS)
{
  dots->push_back(QuantumDot({0.1, 0.1, 0.1}, 0, {0, 0}, {0, 0, 1}));
  dots->push_back(QuantumDot({1.9, 0.1, 0.1}, 0, {0, 0}, {0, 0, 1}));

  AIM::Grid grid(spacing, dots, expansion_order);
  auto expansions =
      Expansions::LeastSquaresExpansionSolver::get_expansions(1, grid, *dots);

  AIM::AimInteraction aim(
      dots, history, interpolation_order, c, dt, grid, expansions,
      AIM::Expansions::Retardation(grid.max_transit_steps(c, dt) +
                                   interpolation_order),
      AIM::normalization::unit);

  for(int t = 0; t < num_steps; ++t) {
    aim.evaluate(t);
  }
}

BOOST_AUTO_TEST_SUITE_END()  // NEARFIELD

BOOST_AUTO_TEST_SUITE_END()  // AIM
