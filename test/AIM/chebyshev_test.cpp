
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>

#include "interactions/AIM/chebyshev.h"
#include "interactions/AIM/projector.h"

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

struct POTENTIAL_PARAMETERS {
  std::vector<QuantumDot> dots;

  Eigen::Vector3d field_fn(const Eigen::Vector3d &r)
  {
    double x = std::pow(r(0), 2) * std::pow(r(1), 2) * std::pow(r(2), 2);
    double y = r.squaredNorm();
    double z = std::cos(2 * M_PI * r(0) / 10) * std::cos(2 * M_PI * r(1) / 10) *
               std::cos(2 * M_PI * r(2) / 10);
    return {x, y, z};
  }
};

BOOST_FIXTURE_TEST_SUITE(POTENTIAL_EVALUATION, POTENTIAL_PARAMETERS)

BOOST_AUTO_TEST_CASE(GRID_SIZE)
{
  constexpr int M = 4;
  std::vector<QuantumDot> dots;
  dots.push_back(QuantumDot({0.1, 0.1, 0.1}));
  dots.push_back(QuantumDot({0.9, 0.9, 0.9}));
  dots.push_back(QuantumDot({1.2, 1.4, 1.8}));

  Eigen::Array<double, 3, 4> spacing;
  spacing.col(0) = Eigen::Vector3d(0.5, 0.5, 0.5);
  spacing.col(1) = Eigen::Vector3d(1, 1, 1);
  spacing.col(2) = Eigen::Vector3d(2, 2, 2);
  spacing.col(3) = Eigen::Vector3d(4, 4, 4);

  auto xs = Chebyshev<M>::roots();
  for(auto &x : xs) x = (x + 1) / 2;

  for(int h = 0; h < 4; ++h) {
    AIM::Grid grid(spacing.col(h), 1, dots);
    std::array<int, 6> shape{{1, grid.size(), M + 1, M + 1, M + 1, 3}};
    boost::multi_array<double, 6> eval(shape), coef(shape);

    std::fill(eval.data(), eval.data() + eval.num_elements(), 0);
    for(int i = 0; i < grid.size(); ++i) {
      const auto r0 = grid.spatial_coord_of_box(i);
      for(int x = 0; x < M + 1; ++x) {
        for(int y = 0; y < M + 1; ++y) {
          for(int z = 0; z < M + 1; ++z) {
            Eigen::Map<Eigen::Vector3d> field(&eval[0][i][x][y][z][0]);
            const Eigen::Vector3d arg =
                r0.array() +
                spacing.col(h) * Eigen::Array3d(xs[x], xs[y], xs[z]);

            field = field_fn(arg);
          }
        }
      }
    }

    Chebyshev<M>::Evaluator<double> bar(grid, dots);
    Chebyshev<M>().fill_coefficients_tensor(grid.size(), eval.data(),
                                            coef.data());

    const boost::multi_array<double, 2> &x =
        bar.interpolate(0, coef, Projector::Potential<double>);

    for(int i = 0; i < static_cast<int>(dots.size()); ++i) {
      Eigen::Map<const Eigen::Vector3d> interp_field(&x[i][0]);
      BOOST_CHECK_SMALL((field_fn(dots[i].position()) - interp_field).norm(),
                        1.2e-3);
      // In general can do much better than this threshold; this just
      // accommodates the large-grid "worst case scenario"
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()  // POTENTIAL_EVALUATION

BOOST_AUTO_TEST_SUITE_END()  // CHEBYSHEV
