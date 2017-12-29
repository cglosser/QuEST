#include "../src/lagrange_set.h"
#include <boost/test/unit_test.hpp>
#include <iostream>

BOOST_AUTO_TEST_SUITE(lagrange_interpolation)

BOOST_AUTO_TEST_CASE(shape_constructor)
{
  constexpr int order = 5;
  Interpolation::UniformLagrangeSet ULS_easy(order);

  BOOST_CHECK(ULS_easy.evaluations.shape()[0] == Interpolation::NUM_DERIVATIVES);
  BOOST_CHECK(ULS_easy.evaluations.shape()[1] == order + 1);
}

BOOST_AUTO_TEST_CASE(value_constructor)
{
  constexpr int order = 5;
  constexpr double eval_point = 0.5;
  Interpolation::UniformLagrangeSet ULS_full(eval_point, order);

  BOOST_CHECK(ULS_full.evaluations.shape()[0] == Interpolation::NUM_DERIVATIVES);
  BOOST_CHECK(ULS_full.evaluations.shape()[1] == order + 1);
}

BOOST_AUTO_TEST_SUITE(mid_point_value_comparison)

BOOST_AUTO_TEST_CASE(full_constructor_d0)
{
  constexpr int order = 5;
  constexpr double eval_point = 0.5;
  const std::vector<double> compare_array = {63.0 / 256,   315.0 / 256,
                                             -105.0 / 128, 63.0 / 128,
                                             -45.0 / 256,  7.0 / 256};

  Interpolation::UniformLagrangeSet ULS_full(eval_point, order);

  BOOST_CHECK_CLOSE(ULS_full.evaluations[0][0], compare_array.at(0), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[0][1], compare_array.at(1), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[0][2], compare_array.at(2), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[0][3], compare_array.at(3), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[0][4], compare_array.at(4), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[0][5], compare_array.at(5), 0.0001);
}

BOOST_AUTO_TEST_CASE(full_constructor_d1)
{
  constexpr int order = 5;
  constexpr double eval_point = 0.5;
  const std::vector<double> compare_array = {-563.0 / 640, 67.0 / 128,
                                             143.0 / 192,  -37.0 / 64,
                                             29.0 / 128,   -71.0 / 1920};

  Interpolation::UniformLagrangeSet ULS_full(eval_point, order);

  BOOST_CHECK_CLOSE(ULS_full.evaluations[1][0], -compare_array.at(0), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[1][1], -compare_array.at(1), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[1][2], -compare_array.at(2), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[1][3], -compare_array.at(3), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[1][4], -compare_array.at(4), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[1][5], -compare_array.at(5), 0.0001);
}

BOOST_AUTO_TEST_CASE(full_constructor_d2)
{
  constexpr int order = 5;
  constexpr double eval_point = 0.5;
  const std::vector<double> compare_array = {
      95.0 / 48, -269.0 / 48, 49.0 / 8, -85.0 / 24, 59.0 / 48, -3.0 / 16};

  Interpolation::UniformLagrangeSet ULS_full(eval_point, order);

  BOOST_CHECK_CLOSE(ULS_full.evaluations[2][0], compare_array.at(0), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[2][1], compare_array.at(1), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[2][2], compare_array.at(2), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[2][3], compare_array.at(3), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[2][4], compare_array.at(4), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[2][5], compare_array.at(5), 0.0001);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(epsilon_value_comparison)

BOOST_AUTO_TEST_CASE(full_constructor_d0)
{
  constexpr int order = 5;
  constexpr double eval_point = 1e-6;
  const std::vector<double> compare_array = {
      39999908666741666638333338333333.0 / 40000000000000000000000000000000.0,
      39999948666690333328666667.0 / 8000000000000000000000000000000.0,
      -19999964333352999995666667.0 / 4000000000000000000000000000000.0,
      13333307333349666662666667.0 / 4000000000000000000000000000000.0,
      -9999979666680333329666667.0 / 8000000000000000000000000000000.0,
      7999983333344999996666667.0 / 40000000000000000000000000000000.0};

  Interpolation::UniformLagrangeSet ULS_full(eval_point, order);

  BOOST_CHECK_CLOSE(ULS_full.evaluations[0][0], compare_array.at(0), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[0][1], compare_array.at(1), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[0][2], compare_array.at(2), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[0][3], compare_array.at(3), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[0][4], compare_array.at(4), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[0][5], compare_array.at(5), 0.0001);
}

BOOST_AUTO_TEST_CASE(full_constructor_d1)
{
  constexpr int order = 5;
  constexpr double eval_point = 1e-6;
  const std::vector<double> compare_array = {
      -18266636666683666662666667.0 / 8000000000000000000000000.0,
      23999938400042599988800001.0 / 4800000000000000000000000.0,
      -3999985733345133329866667.0 / 800000000000000000000000.0,
      2666656266676466663466667.0 / 800000000000000000000000.0,
      -5999975600024599991200001.0 / 4800000000000000000000000.0,
      1599993333340333330666667.0 / 8000000000000000000000000.0};

  Interpolation::UniformLagrangeSet ULS_full(eval_point, order);

  BOOST_CHECK_CLOSE(ULS_full.evaluations[1][0], -compare_array.at(0), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[1][1], -compare_array.at(1), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[1][2], -compare_array.at(2), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[1][3], -compare_array.at(3), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[1][4], -compare_array.at(4), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[1][5], -compare_array.at(5), 0.0001);
}

BOOST_AUTO_TEST_CASE(full_constructor_d2)
{
  constexpr int order = 5;
  constexpr double eval_point = 1e-6;
  const std::vector<double> compare_array = {
      22499974500008999999.0 / 6000000000000000000.0,
      -5133326233336133333.0 / 400000000000000000.0,
      10699982300007799999.0 / 600000000000000000.0,
      -7799985300007199999.0 / 600000000000000000.0,
      2033329233335533333.0 / 400000000000000000.0,
      -4999989500005999999.0 / 6000000000000000000.0};

  Interpolation::UniformLagrangeSet ULS_full(eval_point, order);

  BOOST_CHECK_CLOSE(ULS_full.evaluations[2][0], compare_array.at(0), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[2][1], compare_array.at(1), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[2][2], compare_array.at(2), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[2][3], compare_array.at(3), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[2][4], compare_array.at(4), 0.0001);
  BOOST_CHECK_CLOSE(ULS_full.evaluations[2][5], compare_array.at(5), 0.0001);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
