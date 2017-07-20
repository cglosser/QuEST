#include "../src/integrator/history.h"
#include <Eigen/Dense>
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <memory>
#include "../src/interactions/green_function.h"
#include "../src/interactions/history_interaction.h"
#include "../src/math_utils.h"

BOOST_AUTO_TEST_SUITE(history_interaction)

typedef Eigen::Vector3d vec3d;

struct Universe {
  double e0, c, hbar, dt;
  std::shared_ptr<Propagation::FixedFramePropagator> propagator;

  Universe()
      : e0(3),
        c(2),
        hbar(4),
        dt(0.01),
        propagator(
            std::make_shared<Propagation::FixedFramePropagator>(e0, c, hbar)){};

  Eigen::Vector3d source(double t)
  {
    return Eigen::Vector3d(0, exp(-std::pow(t - 5, 2) / 2.0), 0);
  }

  Eigen::Vector3d mag_d0_source(double t, double delay)
  {
    return Eigen::Vector3d(0, exp(-std::pow(t - 5 - delay, 2) / 2.0), 0);
  }

  Eigen::Vector3d mag_d1_source(double t, double delay)
  {
    return Eigen::Vector3d(
        0, -(t - 5 - delay) * exp(-std::pow(t - 5 - delay, 2) / 2.0), 0);
  }

  Eigen::Vector3d mag_d2_source(double t, double delay)
  {
    return Eigen::Vector3d(0, exp(-std::pow(t - 5 - delay, 2) / 2.0) *
                                  (std::pow(5 + delay - t, 2) - 1),
                           0);
  }

  Eigen::Vector3d analytic_interaction(Eigen::Vector3d &magd1,
                                       Eigen::Vector3d &magd2,
                                       Eigen::Vector3d &dr,
                                       double c,
                                       double e0,
                                       double dist)
  {
    return e0 / (4 * M_PI) *
           dr.cross(
               (magd1 / std::pow(dist, 3) - magd2 / (c * std::pow(dist, 2))));
  }
};

BOOST_FIXTURE_TEST_CASE(history_interaction, Universe)
{
  vec3d pos1(0, 0, 0);
  vec3d pos2(0, 0, 6.5 * c * dt);
  const double total_t = 10;
  const int steps = total_t / dt;
  const double dist = (pos2 - pos1).norm();
  const double delay = dist / c;
  Eigen::Vector3d dr(pos1 - pos2);  // corresponds to separation calculation

  // Set up history with one source "column"
  auto history = std::make_shared<Integrator::History<vec3d>>(2, 22, steps);
  history->fill(Eigen::Vector3d::Zero());

  for(int i = -22; i < steps; ++i) {
    history->array[1][i][0] = source(i * dt);
  }

  // Set up particle list -- don't really care about their "initial" condition
  // (the Eigen:: vector)
  std::shared_ptr<DotVector> dots(std::make_shared<DotVector>(
      DotVector({MagneticParticle(pos1, 1, 1, 1, Eigen::Vector3d::Zero()),
                 MagneticParticle(pos2, 1, 1, 1, Eigen::Vector3d::Zero())})));

  HistoryInteraction history_interaction(dots, history, propagator, 7, dt, c);

  std::cout << std::scientific << std::setprecision(8);
  for(int i = steps * 0.1; i < steps * 0.9; ++i) {
    Eigen::Vector3d magd1 = mag_d1_source(i * dt, delay);
    Eigen::Vector3d magd2 = mag_d2_source(i * dt, delay);
    Eigen::Vector3d interaction =
        analytic_interaction(magd1, magd2, dr, c, e0, dist);

    BOOST_CHECK_CLOSE(interaction[0], history_interaction.evaluate(i)[0][0],
                      1e-6);
  }
}
BOOST_AUTO_TEST_SUITE_END()
