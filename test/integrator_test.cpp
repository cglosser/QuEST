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
  using namespace Integrator::history_enums;

  Integrator::History<int> hist(num_particles, window, num_timesteps);

  BOOST_CHECK_EQUAL(hist.array_.shape()[PARTICLES], num_particles);
  BOOST_CHECK_EQUAL(hist.array_.index_bases()[PARTICLES], 0);

  BOOST_CHECK_EQUAL(hist.array_.shape()[TIMES], num_timesteps + window);
  BOOST_CHECK_EQUAL(hist.array_.index_bases()[TIMES], -window);

  BOOST_CHECK_EQUAL(hist.array_.shape()[DERIVATIVES], 2);
  BOOST_CHECK_EQUAL(hist.array_.index_bases()[DERIVATIVES], 0);

  BOOST_CHECK_EQUAL(hist.array_.num_elements(),
                    num_particles * (num_timesteps + window) * 2);
}

BOOST_AUTO_TEST_CASE(FILLING)
{
  const int fill_value = 4;
  Integrator::History<int> hist(num_particles, window, num_timesteps);
  hist.fill(fill_value);

  for(auto i = 0u; i < hist.array_.num_elements(); ++i) {
    BOOST_CHECK_EQUAL(hist.array_.data()[i], fill_value);
  }
}

BOOST_AUTO_TEST_SUITE_END()

struct SIGMOIDAL_SYSTEM {
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
  double max_abs_error()
  {
    using namespace Integrator::history_enums;

    double max_error = 0;
    for(int i = 1; i < num_steps; ++i) {
      double relative_error =
          (solution(i * dt) - hist->array_[0][i][DERIV_0]) / solution(i * dt);
      max_error = std::max(max_error, std::abs(relative_error));
    }

    return max_error;
  }

  const double dt;
  const int window, num_steps;
  std::shared_ptr<Integrator::History<double>> hist;

  SIGMOIDAL_SYSTEM()
      : dt(0.1),
        window(22),
        num_steps(201),
        hist(
            std::make_shared<Integrator::History<double>>(1, window, num_steps))
  {
    using namespace Integrator::history_enums;

    hist->fill(0);
    for(int i = -window; i <= 0; ++i) {
      hist->array_[0][i][DERIV_0] = solution(i * dt);
      hist->array_[0][i][DERIV_1] = rhs(hist->array_[0][i][DERIV_0], i * dt);
    }
  };
};

BOOST_FIXTURE_TEST_SUITE(ODE_ERROR, SIGMOIDAL_SYSTEM)

BOOST_AUTO_TEST_CASE(PREDICTOR_CORRECTOR)
{
  using namespace Integrator::history_enums;

  std::vector<std::function<double(double, double)>> rhs_funcs{rhs};
  std::unique_ptr<Integrator::RHS<double>> system_rhs =
      std::make_unique<Integrator::ODE_RHS>(dt, hist, rhs_funcs);

  Integrator::PredictorCorrector<double> solver(dt, 32, window, 3.15, hist,
                                                std::move(system_rhs));
  solver.solve();

  BOOST_TEST_MESSAGE(
      "Maximum relative predictor/corrector error: " << max_abs_error());
  BOOST_CHECK_CLOSE(hist->array_[0][num_steps - 1][0],
                    solution((num_steps - 1) * dt), 1e-10);
}

BOOST_AUTO_TEST_CASE(RUNGE_KUTTA_4)
{
  using namespace Integrator::history_enums;

  for(int i = 0; i < num_steps - 1; ++i) {
    const double &yi = hist->array_[0][i][DERIV_0];

    double k1 = rhs(yi, i * dt);
    double k2 = rhs(yi + dt * k1 / 2, (i + 0.5) * dt);
    double k3 = rhs(yi + dt * k2 / 2, (i + 0.5) * dt);
    double k4 = rhs(yi + dt * k3, (i + 1) * dt);

    hist->array_[0][i + 1][DERIV_0] = yi + dt * (k1 + 2 * (k2 + k3) + k4) / 6;
    hist->array_[0][i][DERIV_1] = k1;
  }

  BOOST_TEST_MESSAGE("Maximum relative RK4 error: " << max_abs_error());
  BOOST_CHECK_CLOSE(hist->array_[0][num_steps - 1][0],
                    solution((num_steps - 1) * dt), 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()  // ODE_ERROR

BOOST_AUTO_TEST_SUITE_END()
