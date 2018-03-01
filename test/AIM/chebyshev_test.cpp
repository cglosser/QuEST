
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>

#include "interactions/AIM/chebyshev.h"

BOOST_AUTO_TEST_SUITE(CHEBYSHEV)

BOOST_AUTO_TEST_CASE(QUADRATIC_FIT)
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

BOOST_AUTO_TEST_CASE(FUNCTION_EVALUATION)
{
  std::vector<QuantumDot> dots;
  dots.push_back(QuantumDot({0.1, 0.1, 0.1}));
  dots.push_back(QuantumDot({0.7, 0.8, 0.9}));
  dots.push_back(QuantumDot({1.2, 1.3, 1.4}));
  dots.push_back(QuantumDot({6.2, 7.3, 8.4}));
  AIM::Grid grid(Eigen::Array3d(1, 1, 1), 1, dots);

  constexpr int M = 6;

  std::array<int, 6> shape{{1, grid.size(), M + 1, M + 1, M + 1, 3}};
  boost::multi_array<double, 6> eval(shape), coef(shape);

  const auto field_fn = [](const Eigen::Vector3d &r) -> Eigen::Vector3d {
    double x = std::pow(r(0), 2) * std::pow(r(1), 2) * std::pow(r(2), 2);
    double y = r.squaredNorm();
    double z = std::cos(2 * M_PI * r(0) / 10) * std::cos(2 * M_PI * r(1) / 10) *
               std::cos(2 * M_PI * r(2) / 10);
    return Eigen::Vector3d(x, y, z);
  };

  Chebyshev<M> foo;
  auto xs = foo.roots();

  std::fill(eval.data(), eval.data() + eval.num_elements(), 0);
  for(int i = 0; i < grid.size(); ++i) {
    const auto r0 = grid.spatial_coord_of_box(i);
    for(int x = 0; x < M + 1; ++x) {
      for(int y = 0; y < M + 1; ++y) {
        for(int z = 0; z < M + 1; ++z) {
          const Eigen::Vector3d arg =
              r0.array() +
              spacing * Eigen::Array3d((xs[x] + 1) / 2, (xs[y] + 1) / 2,
                                       (xs[z] + 1) / 2);

          auto f = field_fn(arg);
          eval[0][i][x][y][z][0] = f[0];
          eval[0][i][x][y][z][1] = f[1];
          eval[0][i][x][y][z][2] = f[2];
        }
      }
    }
  }

  foo.fill_coefficients_tensor(grid.size(), eval.data(), coef.data());

  Chebyshev<M>::Evaluator<double> bar(grid, dots);
  const boost::multi_array<double, 2> &x = bar.interpolate(0, coef);

  std::cout.precision(14);
  std::cout << std::scientific << std::endl;
  for(int i = 0; i < static_cast<int>(dots.size()); ++i) {
    std::cout << Eigen::Map<const Eigen::Vector3d>(&x[i][0]).transpose()
              //<< "    " << field_fn(dots[i].position()).transpose()
              << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()  // CHEBYSHEV
