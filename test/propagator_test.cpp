#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../src/interactions/green_function.h"
#include "../src/interactions/rotating_green_function.h"

BOOST_AUTO_TEST_SUITE(propagator)

BOOST_AUTO_TEST_CASE(fixed_frame)
{
  GreenFunction::Dyadic dyadic(4 * M_PI, 1, 1);
  Propagation::FixedFramePropagator ffp(4 * M_PI, 1, 1);

  Eigen::Vector3d dr(1, 1, 1);
  Interpolation::UniformLagrangeSet uls(0.5, 4);

  auto dyadic_mats = dyadic.coefficients(dr, uls);
  auto crtp_mats = ffp.coefficients(dr, uls);

  for(int i = 0; i < 4; ++i) {
    auto delta = dyadic_mats[i] - crtp_mats[i];
    BOOST_CHECK_CLOSE(delta.squaredNorm(), 0, 1e-15);
  }
}

BOOST_AUTO_TEST_CASE(rotating_frame)
{
  GreenFunction::RotatingDyadic rotating_dyadic(4 * M_PI, 1, 1, 2.2);
  Propagation::RotatingFramePropagator rfp(4 * M_PI, 1, 1, 2.2);

  Eigen::Vector3d dr(1, 1, 1);
  Interpolation::UniformLagrangeSet uls(0.5, 4);

  auto rotating_dyadic_mats = rotating_dyadic.coefficients(dr, uls);
  auto rotating_crtp_mats = rfp.coefficients(dr, uls);

  for(int i = 0; i < 4; ++i) {
    auto delta = rotating_dyadic_mats[i] - rotating_crtp_mats[i];
    BOOST_CHECK_CLOSE(delta.squaredNorm(), 0, 1e-15);
  }
}

BOOST_AUTO_TEST_SUITE_END()
