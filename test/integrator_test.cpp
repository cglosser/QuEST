#include "../src/integrator/integrator.h"
#include <boost/test/unit_test.hpp>

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

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
