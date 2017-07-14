#include <Eigen/Dense>
#include <boost/test/unit_test.hpp>
#include <memory>
#include "../src/integrator/history.h"
#include "../src/interactions/green_function.h"
#include "../src/interactions/history_interaction.h"
#include "../src/math_utils.h"

BOOST_AUTO_TEST_SUITE(history_interaction)

typedef Eigen::Vector3d vec3d;

BOOST_AUTO_TEST_CASE(history_interaction)
{
  const double mu0 = 1;
  const double c = 20;
  const double hbar = 1;
  vec3d pos1(0, 0, 0);
  vec3d pos2(0, 0, 0.9);
  double total_t = 10;
  double dt = 0.5e-1;
  int steps = total_t / dt;
  std::vector<Eigen::Vector3d> mag1(steps, vec3d(0, 0, 0)),
      mag2(steps, vec3d(0, 0, 0));
  Eigen::VectorXd times = Eigen::VectorXd::LinSpaced(steps, 0, total_t - dt);

  for(int i = 0; i < steps; ++i) {
    mag1[i][1] = gaussian(times[i] - 5);
    mag2[i][1] = gaussian(times[i] - 4);
  }

  const std::vector<vec3d> mag_control1(mag1);
  const std::vector<vec3d> mag_control2(mag2);

  auto mp1 = MagneticParticle(pos1, 1, 1, 1, mag1[0]);
  auto mp2 = MagneticParticle(pos2, 1, 1, 1, mag2[0]);

  DotVector mpvec = {mp1, mp2};
  auto mpvec_ptr = std::make_shared<DotVector>(mpvec);

  auto dyadic =
      std::make_shared<Propagation::FixedFramePropagator>(mu0, c, hbar);
  
  auto history = std::make_shared<Integrator::History<vec3d>>(2, -22, steps);
 
 
  for(int step = -22; step < steps; ++step) {
    if(step < 0) {
      history->array[0][step][0] = vec3d(0, 0, 0);
      history->array[1][step][0] = vec3d(0, 0, 0);
    } else {
      history->array[0][step][0] = mag1[step];
      history->array[1][step][0] = mag2[step];
    }
  }

  for(int i=-22; i<steps; ++i) std::cout << history->array[0][110][0].transpose() << std::endl;
  
  for(int i=0; i<steps; ++i) {
    BOOST_CHECK(mag1[i] == mag_control1[i]);
    BOOST_CHECK(mag2[i] == mag_control2[i]);
  }

  //for(int i=0; i<steps; ++i) std::cout << mag_control2[i].transpose() << std::endl;
  
  auto hist_inter = HistoryInteraction(mpvec_ptr, history, dyadic, 5, dt, c);
  
  //std::cout << hist_inter.coefficients[0][0] << std::endl;
  
  //for(int i=0; i<steps; ++i) std::cout << hist_inter.evaluate(i)[0].transpose() << std::endl;
}
BOOST_AUTO_TEST_SUITE_END()
