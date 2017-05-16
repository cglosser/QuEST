#include "../src/history.h"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(history)

BOOST_AUTO_TEST_SUITE(nonfinite_detection)

constexpr double dbl_inf = std::numeric_limits<double>::infinity(),
                 dbl_nan = std::numeric_limits<double>::quiet_NaN();
constexpr std::complex<double> zero(0, 0), real_inf(dbl_inf, 0),
    imag_inf(dbl_inf, 0), real_nan(dbl_nan, 0), imag_nan(0, dbl_nan);

BOOST_AUTO_TEST_CASE(zero_is_finite)
{
  History::soltype two_zeros(zero, zero);
  BOOST_CHECK(History::isfinite(two_zeros));
}

BOOST_AUTO_TEST_CASE(real_inf_not_finite)
{
  History::soltype first_is_real_inf(real_inf, zero);
  BOOST_CHECK(!History::isfinite(first_is_real_inf));

  History::soltype second_is_real_inf(zero, real_inf);
  BOOST_CHECK(!History::isfinite(second_is_real_inf));
}

BOOST_AUTO_TEST_CASE(imag_inf_not_finite)
{
  History::soltype first_is_imag_inf(imag_inf, zero);
  BOOST_CHECK(!History::isfinite(first_is_imag_inf));

  History::soltype second_is_imag_inf(zero, imag_inf);
  BOOST_CHECK(!History::isfinite(second_is_imag_inf));
}

BOOST_AUTO_TEST_CASE(real_nan_not_finite)
{
  History::soltype first_is_real_nan(real_nan, zero);
  BOOST_CHECK(!History::isfinite(first_is_real_nan));

  History::soltype second_is_real_nan(zero, real_nan);
  BOOST_CHECK(!History::isfinite(second_is_real_nan));
}

BOOST_AUTO_TEST_CASE(imag_nan_not_finite)
{
  History::soltype first_is_imag_nan(imag_nan, zero);
  BOOST_CHECK(!History::isfinite(first_is_imag_nan));

  History::soltype second_is_imag_nan(zero, imag_nan);
  BOOST_CHECK(!History::isfinite(second_is_imag_nan));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
