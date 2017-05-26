#include "../src/pulse.h"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(pulse_interaction)

BOOST_AUTO_TEST_CASE(pulse_shape_1)
{
  const double amplitude = 15589.2260227;
  const double delay = 5;
  const double width = 227.89013;
  const double freq = 2278.9013;
  const Eigen::Vector3d wavevector(1, 0, 0);
  const Eigen::Vector3d polarization(1, 0, 0);

  const Eigen::Vector3d compare_array(7.7945379681e3,0,0);

  auto pulse = Pulse(amplitude, delay, width, freq, wavevector, polarization);
  auto pulse_eval = pulse(Eigen::Vector3d(1,1,1), 5); 

  BOOST_CHECK_CLOSE(pulse_eval(0), compare_array(0), 1e-6);
  BOOST_CHECK_CLOSE(pulse_eval(1), compare_array(1), 1e-6);
  BOOST_CHECK_CLOSE(pulse_eval(2), compare_array(2), 1e-6);
}

BOOST_AUTO_TEST_CASE(pulse_shape_2)
{
  const double amplitude = 15589.2260227;
  const double delay = 5;
  const double width = 227.89013;
  const double freq = 2278.9013;
  const Eigen::Vector3d wavevector(1, 0, 0);
  const Eigen::Vector3d polarization(1, 0, 0);

  const Eigen::Vector3d compare_array(1.0456587493e3,0,0);

  auto pulse = Pulse(amplitude, delay, width, freq, wavevector, polarization);
  auto pulse_eval = pulse(Eigen::Vector3d(-1,0,1), 5.2); 

  BOOST_CHECK_CLOSE(pulse_eval(0), compare_array(0), 1e-6);
  BOOST_CHECK_CLOSE(pulse_eval(1), compare_array(1), 1e-6);
  BOOST_CHECK_CLOSE(pulse_eval(2), compare_array(2), 1e-6);
}

BOOST_AUTO_TEST_SUITE_END()
