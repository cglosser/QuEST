#include <boost/test/unit_test.hpp>
#include "interactions/AIM/farfield.h"
#include "interactions/AIM/nearfield.h"

BOOST_AUTO_TEST_SUITE(refactor)

struct PARAMETERS {
  int interpolation_order, expansion_order, num_steps, num_dots;
  double c, dt, total_time;

  Eigen::Array3d spacing;

  std::shared_ptr<DotVector> dots;
  std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history;

  PARAMETERS()
      : interpolation_order(3),
        expansion_order(1),
        num_steps(1024),
        num_dots(2),

        c(1),
        dt(1),
        total_time(num_steps * dt),

        spacing(Eigen::Array3d(1, 1, 1) * c * dt),

        dots(std::make_shared<DotVector>()),
        history(std::make_shared<Integrator::History<Eigen::Vector2cd>>(
            num_dots, 10, num_steps))
  {
  }
  double src(const double t) const
  {
    double arg = (t - total_time / 2.0) / (total_time / 12.0);
    return gaussian(arg);
  }

  void fill()
  {
    history->fill(Eigen::Vector2cd::Zero());
    for(int i = -10; i < num_steps; ++i) {
      history->array_[0][i][0](RHO_01) = src(i * dt);
    }
  }
};

BOOST_FIXTURE_TEST_CASE(FARFIELD, PARAMETERS)
{
  dots->push_back(
      QuantumDot(Eigen::Vector3d(0.5, 0.5, 0.5), Eigen::Vector3d(0, 0, 1)));
  dots->push_back(
      QuantumDot(Eigen::Vector3d(0.5, 0.5, 20.5), Eigen::Vector3d(0, 0, 1)));

  fill();

  AIM::Grid grid(spacing, dots, expansion_order);
  const int len = grid.circulant_shape(c, dt, interpolation_order)[0];

  auto expansions =
      AIM::Expansions::LeastSquaresExpansionSolver::get_expansions(
          expansion_order, grid, *dots);

  AIM::Farfield ff(dots, history, interpolation_order, c, dt, grid, expansions,
                   AIM::Expansions::Retardation(len), AIM::normalization::unit);

  std::cout.precision(17);
  for(int t = 0; t < num_steps; ++t) {
    std::cout << t << " " << ff.evaluate(t).transpose() << std::endl;
  }
}

BOOST_FIXTURE_TEST_CASE(NEARFIELD, PARAMETERS)
{
  dots->push_back(
      QuantumDot(Eigen::Vector3d(0.5, 0.5, 0.5), Eigen::Vector3d(0, 0, 1)));
  dots->push_back(
      QuantumDot(Eigen::Vector3d(0.5, 0.5, 1.5), Eigen::Vector3d(0, 0, 1)));

  fill();

  AIM::Grid grid(spacing, dots, expansion_order);
  const int len = grid.circulant_shape(c, dt, interpolation_order)[0];

  auto expansions =
      AIM::Expansions::LeastSquaresExpansionSolver::get_expansions(
          expansion_order, grid, *dots);

  AIM::Nearfield nf(dots, history, interpolation_order, 1, c, dt, grid,
                    expansions, AIM::Expansions::Retardation(len),
                    AIM::normalization::unit);

  AIM::Farfield ff(dots, history, interpolation_order, c, dt, grid,
                    expansions, AIM::Expansions::Retardation(len),
                    AIM::normalization::unit);

  std::cout.precision(17);
  for(int t = 0; t < num_steps; ++t) {
    std::cout << t << " " << nf.evaluate(t).transpose() - ff.evaluate(t).transpose() << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
