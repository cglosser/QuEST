#include "../src/interactions/pulse_interaction.h"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(pulse_interaction)

// Pulse configuration information:
const double amplitude = 15589.2260227;
const double delay = 5;
const double width = 227.89013;
const double pulse_freq = 2278.9013;
const Eigen::Vector3d wavevector(1, 0, 0);
const Eigen::Vector3d polarization(1, 0, 0);

auto pulse =
    Pulse(amplitude, delay, width, pulse_freq, wavevector, polarization);

BOOST_AUTO_TEST_CASE(pulse_shape_1)
{
  const Eigen::Vector3d compare_array(7.7945379681e3, 0, 0);

  auto pulse_eval = pulse(Eigen::Vector3d(1, 1, 1), 5);

  BOOST_CHECK_CLOSE(pulse_eval(0), compare_array(0), 1e-6);
  BOOST_CHECK_CLOSE(pulse_eval(1), compare_array(1), 1e-6);
  BOOST_CHECK_CLOSE(pulse_eval(2), compare_array(2), 1e-6);
}

BOOST_AUTO_TEST_CASE(pulse_shape_2)
{
  const Eigen::Vector3d compare_array(1.0456587493e3, 0, 0);

  auto pulse_eval = pulse(Eigen::Vector3d(-1, 0, 1), 5.2);

  BOOST_CHECK_CLOSE(pulse_eval(0), compare_array(0), 1e-6);
  BOOST_CHECK_CLOSE(pulse_eval(1), compare_array(1), 1e-6);
  BOOST_CHECK_CLOSE(pulse_eval(2), compare_array(2), 1e-6);
}

// QD Configuration Information:
const Eigen::Vector3d pos(1, 1, 1);
const double dot_freq = 2278.9013;
const std::pair<double, double> damping(1, 1);
const Eigen::Vector3d dip(1, 2, 3);

BOOST_AUTO_TEST_CASE(dot_pulse_interaction)
{
  const double compare_value = 1.0641745059e3;

  DotVector dots_vec = {QuantumDot(pos, dot_freq, damping, dip)};
  auto dots = std::make_shared<DotVector>(dots_vec);

  std::shared_ptr<Pulse> pulse_ptr = std::make_shared<Pulse>(
      Pulse(amplitude, delay, width, pulse_freq, wavevector, polarization));

  PulseInteraction pulse_interaction =
      PulseInteraction(dots, pulse_ptr, 1, 0.1);  // hbar=1, dt=1

  auto results = pulse_interaction.evaluate(52);

  BOOST_CHECK_CLOSE(real(results(0)), compare_value, 1e-6);
}
BOOST_AUTO_TEST_SUITE_END()
