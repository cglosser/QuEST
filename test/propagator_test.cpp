#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../src/interactions/green_function.h"

BOOST_AUTO_TEST_SUITE(propagator)

BOOST_AUTO_TEST_CASE(fixed_frame)
{
  Propagation::FixedFramePropagator ffp(4 * M_PI, 1, 1);

  Interpolation::UniformLagrangeSet uls(0.5, 4);
  auto mats = ffp.coefficients(Eigen::Vector3d(1, 1, 1), uls);

  for(const auto &m : mats) {
    std::cout << m << std::endl << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(rotating_frame)
{
  Propagation::RotatingFramePropagator rfp(4 * M_PI, 1, 1, 2.2);

  Interpolation::UniformLagrangeSet uls(0.5, 4);
  auto mats = rfp.coefficients(Eigen::Vector3d(1, 1, 1), uls);

  for(const auto &m : mats) {
    std::cout << m << std::endl << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
