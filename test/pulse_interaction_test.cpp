#include "../src/interactions/pulse_interaction.h"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(pulse_interaction)

struct Universe {
  double amp, delay, width, freq, dt;
  Eigen::Vector3d wavevector, magnetization;
  Pulse pulse;

  Universe()
      : amp(2),
        delay(0.5),
        width(2),
        freq(10),
        dt(1),
        wavevector(1, 1, 1),
        magnetization(1, 1, 1),
        pulse(Pulse(amp, delay, width, freq, wavevector, magnetization)){};

  Eigen::Vector3d evaluate(double amp,
                           double delay,
                           double width,
                           double freq,
                           double t,
                           Eigen::Vector3d wavevec,
                           Eigen::Vector3d mag,
                           Eigen::Vector3d pos)
  {
    double arg = wavevec.dot(pos) - freq * (t - delay);
    return (amp * 2 * mag.normalized() *
            std::exp(-std::pow(arg / width, 2) / 2) * cos(arg));
  }
};

BOOST_FIXTURE_TEST_CASE(pulse_shape, Universe)
{
  const double t = 1;
  const Eigen::Vector3d pos(1, 1, 1);

  const Eigen::Vector3d compare_array(
      evaluate(amp, delay, width, freq, t, wavevector, magnetization, pos));

  auto pulse_eval = pulse(pos, t);

  BOOST_CHECK_CLOSE(pulse_eval(0), compare_array(0), 1e-15);
  BOOST_CHECK_CLOSE(pulse_eval(1), compare_array(1), 1e-15);
  BOOST_CHECK_CLOSE(pulse_eval(2), compare_array(2), 1e-15);
}

BOOST_FIXTURE_TEST_CASE(dot_pulse_interaction, Universe)
{
  const Eigen::Vector3d pos(1, 1, 1);
  const Eigen::Vector3d particle_mag(1, 0, 1);
  const double alpha = 2;
  const double gamma0 = 2;
  const double sat_mag = 2;
  const double t = 1;

  DotVector dots_vec = {
      MagneticParticle(pos, alpha, gamma0, sat_mag, particle_mag)};
  auto dots = std::make_shared<DotVector>(dots_vec);

  std::shared_ptr<Pulse> pulse_ptr = std::make_shared<Pulse>(pulse);

  PulseInteraction pulse_interaction = PulseInteraction(dots, pulse_ptr, 1, dt);

  auto results = pulse_interaction.evaluate(1);
  auto analytic_pulse =
      evaluate(amp, delay, width, freq, t, wavevector, magnetization, pos);

  BOOST_CHECK(results[0][0] == analytic_pulse[0]);
  BOOST_CHECK(results[0][1] == analytic_pulse[1]);
  BOOST_CHECK(results[0][2] == analytic_pulse[2]);
}
BOOST_AUTO_TEST_SUITE_END()
