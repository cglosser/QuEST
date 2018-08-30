#include <boost/test/unit_test.hpp>

#include <rapidcheck/boost_test.h>
#include <iostream>
#include "../src/math_utils.h"

BOOST_AUTO_TEST_SUITE(math_utils)

BOOST_AUTO_TEST_SUITE(LINSPACE)

RC_BOOST_PROP(first_arg_is_first_element,
              (const double first, const double second))
{
  const auto values = linspace(first, second, 8);
  const auto error = std::abs(values.at(0) - first);
  RC_ASSERT(error <= 1e-12 * std::abs(first));
  RC_ASSERT(error <= 1e-12 * std::abs(values.at(0)));
}

RC_BOOST_PROP(second_arg_is_last_element,
              (const double first, const double second))
{
  const auto values = linspace(first, second, 8);
  const auto error = std::abs(values.at(7) - second);
  RC_ASSERT(error <= 1e-12 * std::abs(second));
  RC_ASSERT(error <= 1e-12 * std::abs(values.at(7)));
}

RC_BOOST_PROP(reverse_args_reverses_list,
              (const double first, const double second))
{
  constexpr size_t n = 8;

  auto l1 = linspace(first, second, n);
  auto l2 = linspace(second, first, n);
  std::reverse(l2.begin(), l2.end());

  for(size_t i = 0; i < 8; ++i) {
    const double error = std::abs(l1.at(i) - l2.at(i));
    RC_ASSERT(error <= 1e-12 * std::abs(l1.at(i)));
    RC_ASSERT(error <= 1e-12 * std::abs(l2.at(i)));
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(UNIT_NORMAL)

RC_BOOST_PROP(vector_has_unit_norm, (const double theta, const double phi))
{
  const Eigen::Vector3d r = unit_normal(theta, phi);
  RC_ASSERT(std::abs(r.norm() - 1) <= 1e-12);
}

RC_BOOST_PROP(vector_has_bounded_coordinates,
              (const double theta, const double phi))
{
  const Eigen::Vector3d r = unit_normal(theta, phi);

  RC_ASSERT(-1 <= r(0) && r(0) <= 1);
  RC_ASSERT(-1 <= r(1) && r(1) <= 1);
  RC_ASSERT(-1 <= r(2) && r(2) <= 1);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GRID_SEQUENCE)

BOOST_AUTO_TEST_CASE(GridSequence)
{
  BOOST_CHECK_EQUAL(grid_sequence(0), 0);
  BOOST_CHECK_EQUAL(grid_sequence(1), 1);
  BOOST_CHECK_EQUAL(grid_sequence(2), -1);
  BOOST_CHECK_EQUAL(grid_sequence(3), 2);
  BOOST_CHECK_EQUAL(grid_sequence(4), -2);
  BOOST_CHECK_EQUAL(grid_sequence(5), 3);
  BOOST_CHECK_EQUAL(grid_sequence(6), -3);
  BOOST_CHECK_EQUAL(grid_sequence(7), 4);
  BOOST_CHECK_EQUAL(grid_sequence(8), -4);
  BOOST_CHECK_EQUAL(grid_sequence(9), 5);
}

RC_BOOST_PROP(only_even_terms_are_negative, ())
{
  const auto n = *rc::gen::inRange(0, 4096);
  const auto x = grid_sequence(n);

  if(n % 2 == 0) {
    RC_ASSERT(x <= 0);
  } else {
    RC_ASSERT(x > 0);
  }
}

RC_BOOST_PROP(grid_sequence_pattern, ())
{
  const auto n = *rc::gen::inRange(0, 4096);
  const auto x = grid_sequence(n);

  if(x < 0) {
    RC_ASSERT(grid_sequence(n + 1) == -x + 1);
  } else if(x > 0) {
    RC_ASSERT(grid_sequence(n + 1) == -x);
  } else {
    RC_ASSERT(grid_sequence(n + 1) == 1);
  }
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SPLIT_DOUBLE)

RC_BOOST_PROP(int_magnitude_is_smaller_than_double, (const double x))
{
  const auto p = split_double(x);
  RC_ASSERT(std::abs(p.first) <= std::abs(x));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(FALLING_FACTORIAL)
{
  BOOST_CHECK_EQUAL(falling_factorial(1, 0), 1);
  BOOST_CHECK_EQUAL(falling_factorial(1.5, 0), 1);
  BOOST_CHECK_EQUAL(falling_factorial(2, 0), 1);
  BOOST_CHECK_EQUAL(falling_factorial(3.14159, 0), 1);

  BOOST_CHECK_EQUAL(falling_factorial(1, 1), 1);
  BOOST_CHECK_EQUAL(falling_factorial(1, 2), 0);
  BOOST_CHECK_EQUAL(falling_factorial(1, 3), 0);

  BOOST_CHECK_EQUAL(falling_factorial(1, 1), 1);
  BOOST_CHECK_EQUAL(falling_factorial(2, 1), 2);
  BOOST_CHECK_EQUAL(falling_factorial(3, 1), 3);

  BOOST_CHECK_EQUAL(falling_factorial(1.5, 2), 3.0 / 4.0);
  BOOST_CHECK_EQUAL(falling_factorial(2.5, 2), 15.0 / 4.0);
  BOOST_CHECK_EQUAL(falling_factorial(3.5, 2), 35.0 / 4.0);
}

BOOST_AUTO_TEST_SUITE_END()
