// clang-format off
#include "../src/history.h"
#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <vector>
// clang-format on

namespace data = boost::unit_test::data;

BOOST_AUTO_TEST_SUITE(history)

BOOST_AUTO_TEST_SUITE(nonfinite_detection)

constexpr double dbl_inf = std::numeric_limits<double>::infinity(),
                 dbl_nan = std::numeric_limits<double>::quiet_NaN();
constexpr std::complex<double> zero(0, 0);

BOOST_AUTO_TEST_CASE(zero_is_finite)
{
  History::soltype two_zeros(zero, zero);
  BOOST_CHECK(History::isfinite(two_zeros));
}

const std::vector<std::complex<double>> nonfinite_values = {
    std::complex<double>(dbl_inf, 0), std::complex<double>(0, dbl_inf),
    std::complex<double>(dbl_nan, 0), std::complex<double>(0, dbl_nan),
};

BOOST_DATA_TEST_CASE(nonfinite_detection, nonfinite_values)
{
  History::soltype first(sample, zero), second(zero, sample);
  BOOST_TEST(!History::isfinite(first));
  BOOST_TEST(!History::isfinite(second));
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
