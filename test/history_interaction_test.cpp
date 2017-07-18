#include <Eigen/Dense>
#include <boost/test/unit_test.hpp>
#include <memory>
#include <fstream>
#include "../src/integrator/history.h"
#include "../src/interactions/green_function.h"
#include "../src/interactions/history_interaction.h"
#include "../src/math_utils.h"

BOOST_AUTO_TEST_SUITE(history_interaction)

typedef Eigen::Vector3d vec3d;

struct Universe {
  double mu0, c, hbar, dt;
  std::shared_ptr<Propagation::FixedFramePropagator> propagator;

  Universe()
      : mu0(1),
        c(1),
        hbar(1),
        dt(0.05),
        propagator(std::make_shared<Propagation::FixedFramePropagator>(
            mu0, c, hbar)){};

  Eigen::Vector3d source(double t)
  {
    return Eigen::Vector3d(0, exp(-std::pow(t - 5, 2) / 2.0), 0);
  }
};

BOOST_FIXTURE_TEST_CASE(history_interaction, Universe)
{
  vec3d pos1(0, 0, 0);
  vec3d pos2(0, 0, 0.5 * c * dt);
  const double total_t = 10;
  const int steps = total_t / dt;

  // Set up history with one source "column"
  auto history = std::make_shared<Integrator::History<vec3d>>(2, 22, steps);
  history->fill(Eigen::Vector3d::Zero());

  for(int i = -22; i < steps; ++i) {
    history->array[1][i][0] = source(i * dt);
  }
  
  // Set up particle list -- don't really care about their "initial" condition (the Eigen:: vector)
  std::shared_ptr<DotVector> dots(std::make_shared<DotVector>(
      DotVector({MagneticParticle(pos1, 1, 1, 1, Eigen::Vector3d::Zero()),
                 MagneticParticle(pos2, 1, 1, 1, Eigen::Vector3d::Zero())})));


  HistoryInteraction history_interaction(dots, history, propagator, 5, dt, c);  

  std::cout << std::scientific << std::setprecision(8);
  for(int i=0; i<steps; ++i) {
    std::cout << i << " ";
    std::cout << history_interaction.evaluate(i)[0].transpose() << " | ";
    std::cout << history_interaction.evaluate(i)[0].transpose() << std::endl;
  }
}
BOOST_AUTO_TEST_SUITE_END()
