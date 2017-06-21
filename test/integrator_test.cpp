#include "../src/integrator/integrator.h"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(integrator)

BOOST_AUTO_TEST_SUITE(history)

BOOST_AUTO_TEST_CASE(templated_types)
{
  Integrator::History<int> int_hist(0, 0, 0);
  Integrator::History<double> dbl_hist(0, 0, 0);
  Integrator::History<Eigen::Vector2cd> cvec_hist(0, 0, 0);
}

BOOST_AUTO_TEST_CASE(shape)
{
  const int num_particles = 2;
  const int window = 4;
  const int num_timesteps = 8;
  Integrator::History<int> hist(num_particles, window, num_timesteps);

  BOOST_CHECK(hist.array.shape()[0] == num_particles);
  BOOST_CHECK(hist.array.shape()[1] == num_timesteps + window);
  BOOST_CHECK(hist.array.shape()[2] == 2);

  BOOST_CHECK(hist.array.index_bases()[0] == 0);
  BOOST_CHECK(hist.array.index_bases()[1] == -window);
  BOOST_CHECK(hist.array.index_bases()[2] == 0);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
