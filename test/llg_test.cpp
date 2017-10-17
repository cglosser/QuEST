#include <Eigen/Dense>
#include <boost/test/unit_test.hpp>
#include "../src/magnetic_particle.h"

BOOST_AUTO_TEST_SUITE(llg_rhs)

BOOST_AUTO_TEST_CASE(llg_rhs)
{
  typedef Eigen::Vector3d vec3d;

  vec3d pos(1, 2, 1);
  const vec3d mag(100, 100, 100);
  const double alpha = 2;
  const double gamma0 = 3;
  const double sat_mag = 100;

  MagneticParticle mp(pos, alpha, gamma0, sat_mag, mag);
  const vec3d field(10,100,50);

  vec3d mp.llg_rhs(mag, field);

}
BOOST_AUTO_TEST_SUITE_END()
