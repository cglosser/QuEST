#include <boost/test/unit_test.hpp>
#include <iomanip>

#include "interactions/AIM/aim_interaction.h"
#include "interactions/AIM/nearfield_interaction.h"
#include "interactions/direct_interaction.h"

BOOST_AUTO_TEST_SUITE(DIRECT_AIM_COMPARISON)

struct PARAMETERS {
  int interpolation_order, expansion_order, num_steps, num_dots,
      prehistory_length;
  double c, dt, total_time, omega;

  Eigen::Array3d spacing;

  std::shared_ptr<DotVector> dots;
  std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history;

  PARAMETERS(const int num_dots)
      : interpolation_order(4),
        expansion_order(4),
        num_steps(1024),
        num_dots(num_dots),
        prehistory_length(10),

        c(1),
        dt(1),
        total_time(dt * num_steps),
        omega(M_PI / 10),

        spacing(Eigen::Array3d(1, 1, 1) * c * dt),

        dots(std::make_shared<DotVector>(num_dots)),
        history(std::make_shared<Integrator::History<Eigen::Vector2cd>>(
            num_dots, prehistory_length, num_steps))
  {
    history->fill(Eigen::Vector2cd::Zero());
    for(int d = 0; d < num_dots; ++d) {
      for(int i = -10; i < num_steps; ++i) {
        // Everybody radiates a Gaussian
        history->array[d][i][0](RHO_01) = src(i * dt);
      }
    }
  };

  double src(const double t) const
  {
    double arg = (t - total_time / 2.0) / (total_time / 6.0);
    return gaussian(arg);
  }
};

struct TWO_PARTICLE : public PARAMETERS {
  TWO_PARTICLE() : PARAMETERS(2){};
};

struct THREE_PARTICLE : public PARAMETERS {
  THREE_PARTICLE() : PARAMETERS(3){};
};

BOOST_FIXTURE_TEST_CASE(NEARFIELD, TWO_PARTICLE)
{
  const Eigen::Vector3d zhat(0, 0, 1);
  dots->at(0) = QuantumDot(Eigen::Vector3d(0.1, 0.1, 0.1), zhat);
  dots->at(1) = QuantumDot(Eigen::Vector3d(0.9, 0.9, 0.9), zhat);

  Propagation::RotatingFramePropagator greens_function(1, c, omega);

  DirectInteraction direct(dots, history, greens_function, interpolation_order,
                           c, dt);

  AIM::Grid grid(spacing, dots, expansion_order);
  AIM::NearfieldInteraction nf(dots, history, greens_function,
                               interpolation_order, c, dt, grid);

  for(int t = 0; t < num_steps; ++t) {
    BOOST_CHECK_EQUAL((direct.evaluate(t) - nf.evaluate(t)).matrix().norm(), 0);
  }
}

BOOST_FIXTURE_TEST_CASE(FARFIELD, TWO_PARTICLE)
{
  const Eigen::Vector3d zhat(0, 0, 1);
  dots->at(0) = QuantumDot(Eigen::Vector3d(0.1, 0.1, 0.1), zhat);
  dots->at(1) = QuantumDot(Eigen::Vector3d(9.9, 9.9, 9.9), zhat);

  Propagation::RotatingFramePropagator greens_function(1, c, omega);

  DirectInteraction direct(dots, history, greens_function, interpolation_order,
                           c, dt);

  AIM::Grid grid(spacing, dots, expansion_order);
  AIM::NearfieldInteraction nf(dots, history, greens_function,
                               interpolation_order, c, dt, grid);

  BOOST_CHECK_EQUAL(nf.build_pair_list().size(), 0);

  auto expansion_table =
      AIM::Expansions::LeastSquaresExpansionSolver::get_expansions(
          expansion_order, grid, *dots);
  AIM::AimInteraction aim(
      dots, history, interpolation_order, c, dt, grid, expansion_table,
      AIM::Expansions::RotatingEFIE(
          grid.max_transit_steps(c, dt) + interpolation_order, c, dt, omega),
      AIM::normalization::Helmholtz(omega / c, 1));

  for(int t = 0; t < num_steps; ++t) {
    BOOST_CHECK_SMALL((direct.evaluate(t) - aim.evaluate(t)).matrix().norm(),
                      1e-2);
  }
}

BOOST_FIXTURE_TEST_CASE(NO_AIM_NEARFIELD, THREE_PARTICLE)
{
  const Eigen::Vector3d zhat(0, 0, 1);
  dots->at(0) = QuantumDot(Eigen::Vector3d(0.1, 0.1, 0.1), zhat);
  dots->at(1) = QuantumDot(Eigen::Vector3d(0.9, 0.9, 0.9), zhat);
  dots->at(2) = QuantumDot(Eigen::Vector3d(9.9, 9.9, 9.9), zhat);

  for(int t = 0; t < num_steps; ++t) {
    history->array[1][t][0][1] = history->array[2][t][0][1] = 0;
  }

  Propagation::RotatingFramePropagator greens_function(1, c, omega);

  DirectInteraction direct(dots, history, greens_function, interpolation_order,
                           c, dt);

  AIM::Grid grid(spacing, dots, expansion_order);
  AIM::NearfieldInteraction nf(dots, history, greens_function,
                               interpolation_order, c, dt, grid);

  BOOST_CHECK_EQUAL(nf.build_pair_list().size(), 1);

  auto expansion_table =
      AIM::Expansions::LeastSquaresExpansionSolver::get_expansions(
          expansion_order, grid, *dots);
  AIM::AimInteraction aim(
      dots, history, interpolation_order, c, dt, grid, expansion_table,
      AIM::Expansions::RotatingEFIE(
          grid.max_transit_steps(c, dt) + interpolation_order, c, dt, omega),
      AIM::normalization::Helmholtz(omega / c, 1));

  for(int t = 0; t < num_steps; ++t) {
    Eigen::ArrayXcd aim_eval = aim.evaluate(t);
    BOOST_CHECK_SMALL(std::abs(aim_eval(1)),
                      1e-15);  // No nearfield component from AIM
    BOOST_CHECK_SMALL(std::abs(direct.evaluate(t)(2) - aim_eval(2)),
                      1e-3);  // Agrees with direct
  }
}

BOOST_AUTO_TEST_SUITE_END()  // DIRECT_AIM_COMPARISON
