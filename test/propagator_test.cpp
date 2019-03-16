#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../src/integrator/history.h"
#include "../src/interactions/direct_interaction.h"
#include "../src/interactions/green_function.h"
#include "../src/math_utils.h"
#include "../src/quantum_dot.h"

BOOST_AUTO_TEST_SUITE(PROPAGATION)

struct Parameters {
  const double c0;
  const double dt;
  const int num_timesteps;
  const int interpolation_order;

  const double max_time;

  Parameters()
      : c0{1},
        dt{1},
        num_timesteps{2048},
        interpolation_order{6},
        max_time{num_timesteps * dt} {};
};

BOOST_AUTO_TEST_SUITE(DIRECT_INTERACTION)

BOOST_FIXTURE_TEST_CASE(RETARDED_GAUSSIAN, Parameters)
{
  // == Initialize dots ==
  auto qds = std::make_shared<DotVector>();
  qds->emplace_back(Eigen::Vector3d(0.5, 0.5, 0.5), Eigen::Vector3d(0, 0, 1));
  qds->emplace_back(Eigen::Vector3d(10.5, 0.6, 0.7), Eigen::Vector3d(0, 0, 1));

  const double tau =
      (qds->at(1).position() - qds->at(0).position()).norm() / c0;

  // == Initialize history ==
  auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(
      2, num_timesteps / 2, num_timesteps);

  const double mu = max_time / 2.0, sigma = max_time / 12.0;

  for(int t = -num_timesteps / 2; t < num_timesteps; ++t) {
    const double arg = (t * dt - mu) / sigma;

    history->array_[0][t][0] = Eigen::Vector2cd(0, gaussian(arg));
    history->array_[1][t][0] = Eigen::Vector2cd(0, 0);
  }

  // == Initialize interaction ==
  Propagation::Identity<cmplx> identity_kernel;
  DirectInteraction direct_interaction(qds, history, identity_kernel,
                                       interpolation_order, c0, dt);

  // == Run propagation simulation ==
  for(int t = 0; t < num_timesteps; ++t) {
    BOOST_TEST_MESSAGE(t);

    const double obs_val_calculated = direct_interaction.evaluate(t)(1).real();
    const double obs_val_actual = gaussian((t * dt - mu - tau) / sigma);

    BOOST_REQUIRE_CLOSE_FRACTION(obs_val_calculated, obs_val_actual, 1e-12);
  }
}

BOOST_AUTO_TEST_SUITE_END()  // DIRECT_INTERACTION

BOOST_AUTO_TEST_SUITE_END()  // PROPAGATION
