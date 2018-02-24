#include <boost/multi_array.hpp>
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>

#include "interactions/AIM/chebyshev.h"

BOOST_AUTO_TEST_SUITE(CHEBYSHEV)

BOOST_AUTO_TEST_CASE(QUADRATIC)
{
  constexpr int num_boxes = 5, M = 2;

  constexpr std::array<int, 5> shape{{num_boxes, M + 1, M + 1, M + 1, 3}};
  boost::multi_array<double, 5> eval(shape), coef(shape);
  std::fill(eval.data(), eval.data() + eval.num_elements(), 0);

  auto xs = Chebyshev<M>::roots();

  for(int x = 0; x < M + 1; ++x) {
    for(int y = 0; y < M + 1; ++y) {
      for(int z = 0; z < M + 1; ++z) {
        eval[0][x][y][z][0] = std::pow(xs[x], 2);
        eval[1][x][y][z][0] = std::pow(xs[y], 2);
        eval[2][x][y][z][0] = std::pow(xs[z], 2);
        eval[3][x][y][z][0] =
            std::pow(xs[x], 2) + std::pow(xs[y], 2) + std::pow(xs[z], 2);
        eval[4][x][y][z][0] =
            std::pow(xs[x], 2) * std::pow(xs[y], 2) * std::pow(xs[z], 2);
      }
    }
  }

  Chebyshev<M> cheb;
  cheb.fill_coefficients_tensor(num_boxes, eval.data(), coef.data());
  for(size_t i = 0; i < coef.num_elements(); ++i) {
    if(std::abs(*(coef.data() + i)) < 1e-12) *(coef.data() + i) = 0;
  }

  for(int x = 0; x < M + 1; ++x) {
    if(x == 0 || x == 2) BOOST_CHECK_CLOSE(coef[0][x][0][0][0], 0.5, 1e-12);
    for(int y = 0; y < M + 1; ++y) {
      if(y == 0 || y == 2) BOOST_CHECK_CLOSE(coef[1][0][y][0][0], 0.5, 1e-12);
      for(int z = 0; z < M + 1; ++z) {
        if(z == 0 || z == 2) BOOST_CHECK_CLOSE(coef[2][0][0][z][0], 0.5, 1e-12);

        BOOST_CHECK_CLOSE(
            coef[3][x][y][z][0],
            coef[0][x][y][z][0] + coef[1][x][y][z][0] + coef[2][x][y][z][0],
            1e-12);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()  // CHEBYSHEV
