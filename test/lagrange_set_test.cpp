#include <boost/test/unit_test.hpp>
#include "../src/lagrange_set.h"

BOOST_AUTO_TEST_SUITE(lagrange_interpolation)

BOOST_AUTO_TEST_CASE(shape_constructor)
{
  constexpr int order = 5;
  Interpolation::UniformLagrangeSet ULS_easy(order);

  BOOST_CHECK(ULS_easy.weights.shape()[0] == Interpolation::NUM_DERIVATIVES);
  BOOST_CHECK(ULS_easy.weights.shape()[1] == order + 1);
}

BOOST_AUTO_TEST_CASE(value_constructor)
{
  constexpr int order = 5;
  constexpr double eval_point = 0.5;
  Interpolation::UniformLagrangeSet ULS_full(eval_point, order);

  BOOST_CHECK(ULS_full.weights.shape()[0] == Interpolation::NUM_DERIVATIVES);
  BOOST_CHECK(ULS_full.weights.shape()[1] == order + 1);
}

BOOST_AUTO_TEST_SUITE_END()
