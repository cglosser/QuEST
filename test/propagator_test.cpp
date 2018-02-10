#include <boost/test/unit_test.hpp>
#include <iostream>

#include "../src/interactions/green_function.h"

BOOST_AUTO_TEST_SUITE(propagator)

BOOST_AUTO_TEST_CASE(fixed_frame)
{
  Propagation::EFIE<double> ffp(1, 1);

  Eigen::Vector3d dr(1, 1, 1);
  Interpolation::UniformLagrangeSet uls(0.5, 4);

  auto crtp_mats = ffp.coefficients(dr, uls);

  //for(const auto &m : crtp_mats) {
    //std::cout << m << std::endl << std::endl;
  //}
}

BOOST_AUTO_TEST_CASE(rotating_frame)
{
  Propagation::RotatingEFIE rfp(1, 1, 2.2);

  Eigen::Vector3d dr(1, 1, 1);
  Interpolation::UniformLagrangeSet uls(0.5, 4);

  auto rotating_crtp_mats = rfp.coefficients(dr, uls);

  //for(const auto &m : rotating_crtp_mats) {
    //std::cout << m << std::endl << std::endl;
  //}
}

BOOST_AUTO_TEST_SUITE_END()
