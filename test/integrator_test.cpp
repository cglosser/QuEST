#include <boost/test/unit_test.hpp>
#include <cmath>

#include "../src/integrator/RHS/ode_rhs.h"
#include "../src/integrator/integrator.h"

BOOST_AUTO_TEST_SUITE(INTEGRATOR)

struct DIMENSIONS {
  int num_particles, window, num_timesteps;
  DIMENSIONS() : num_particles(2), window(4), num_timesteps(8){};
};

BOOST_FIXTURE_TEST_SUITE(HISTORY, DIMENSIONS)

BOOST_AUTO_TEST_CASE(CONSTRUCTION)
{
  using namespace Integrator::enums;

  Integrator::History<int> hist(num_particles, window, num_timesteps);

  BOOST_CHECK_EQUAL(hist.array.shape()[PARTICLES], num_particles);
  BOOST_CHECK_EQUAL(hist.array.index_bases()[PARTICLES], 0);

  BOOST_CHECK_EQUAL(hist.array.shape()[TIMES], num_timesteps + window);
  BOOST_CHECK_EQUAL(hist.array.index_bases()[TIMES], -window);

  BOOST_CHECK_EQUAL(hist.array.shape()[DERIVATIVES], 2);
  BOOST_CHECK_EQUAL(hist.array.index_bases()[DERIVATIVES], 0);

  BOOST_CHECK_EQUAL(hist.array.num_elements(),
                    num_particles * (num_timesteps + window) * 2);
}

BOOST_AUTO_TEST_CASE(FILLING)
{
  const int fill_value = 4;
  Integrator::History<int> hist(num_particles, window, num_timesteps);
  hist.fill(fill_value);

  for(auto i = 0u; i < hist.array.num_elements(); ++i) {
    BOOST_CHECK_EQUAL(hist.array.data()[i], fill_value);
  }
}

BOOST_AUTO_TEST_SUITE_END()

struct SigmoidalSystem {
  static double rhs(double f, double t)
  {
    return 1 / (1 + std::exp(-(t - 10))) - f;
  }
  static double solution(double t)
  {
    return (-1 + std::exp(t) + std::exp(10) * std::log(1 + std::exp(10)) -
            std::exp(10) * std::log(std::exp(10) + std::exp(t))) /
           std::exp(t);
  }
};

BOOST_FIXTURE_TEST_CASE(ODE_ERROR, SigmoidalSystem)
{
  const double dt = 0.1;
  auto hist = std::make_shared<Integrator::History<double>>(1, 22, 201);
  std::vector<std::function<double(double, double)>> rhs_funcs{rhs};
  std::unique_ptr<Integrator::RHS<double>> system_rhs =
      std::make_unique<Integrator::ODE_RHS>(dt, hist, rhs_funcs);

  hist->fill(0);
  for(int i = -22; i <= 0; ++i) {
    hist->array[0][i][0] = solution(i * dt);
    hist->array[0][i][1] = rhs(hist->array[0][i][0], i * dt);
  }

  Integrator::PredictorCorrector<double> solver(dt, 18, 22, 3.15, hist,
                                                std::move(system_rhs));
  solver.solve();

  BOOST_CHECK_CLOSE(hist->array[0][200][0], solution(200 * dt), 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()
