#include <boost/test/unit_test.hpp>
#include <cmath>

#include "../src/integrator/RHS/rhs.h"
#include "../src/integrator/integrator.h"

BOOST_AUTO_TEST_SUITE(integrator)

struct Shape {
  int num_particles, window, num_timesteps;
  Shape() : num_particles(2), window(4), num_timesteps(8){};
};

BOOST_FIXTURE_TEST_SUITE(history, Shape)

BOOST_AUTO_TEST_CASE(templated_types)
{
  Integrator::History<int> int_hist(num_particles, window, num_timesteps);
  Integrator::History<double> dbl_hist(num_particles, window, num_timesteps);
  Integrator::History<Eigen::Vector2cd> cvec_hist(num_particles, window,
                                                  num_timesteps);
}

BOOST_AUTO_TEST_CASE(shape)
{
  Integrator::History<int> hist(num_particles, window, num_timesteps);

  BOOST_CHECK(hist.array.shape()[0] == static_cast<size_t>(num_particles));
  BOOST_CHECK(hist.array.shape()[1] ==
              static_cast<size_t>(num_timesteps + window));
  BOOST_CHECK(hist.array.shape()[2] == 2);

  BOOST_CHECK(hist.array.index_bases()[0] == 0);
  BOOST_CHECK(hist.array.index_bases()[1] == -window);
  BOOST_CHECK(hist.array.index_bases()[2] == 0);
}

BOOST_AUTO_TEST_CASE(filling)
{
  const int fill_value = 4;
  Integrator::History<int> hist(num_particles, window, num_timesteps);
  hist.fill(fill_value);

  for(int n = 0; n < num_particles; ++n) {
    for(int t = -window; t < num_timesteps; ++t) {
      BOOST_CHECK(hist.array[n][t][0] == fill_value);
      BOOST_CHECK(hist.array[n][t][1] == fill_value);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()

struct SigmoidalSystem {
  static double rhs(double f, double t) { return 1 / (1 + std::exp(-(t - 10))) - f; }
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
      std::make_unique<Integrator::RHS<double>>(dt, hist, rhs_funcs);

  hist->fill(0);
  for(int i = -22; i <= 0; ++i) {
    hist->array[0][i][0] = solution(i * dt);
    hist->array[0][i][1] = rhs(hist->array[0][i][0], i * dt);
  }

  Integrator::PredictorCorrector<double> solver(dt, 18, 22, 3.15, hist,
                                                system_rhs);
  solver.solve();

  BOOST_CHECK_CLOSE(hist->array[0][200][0], solution(20), 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()
