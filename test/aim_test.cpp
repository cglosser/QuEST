#include <boost/test/unit_test.hpp>
#include <cmath>

#include "../src/interactions/aim_interaction.h"

BOOST_AUTO_TEST_SUITE(AIM)

BOOST_AUTO_TEST_SUITE(GridTest)

BOOST_AUTO_TEST_SUITE(TwoDimensions)

struct PointSetup {
  std::pair<double, double> damping;
  Eigen::Vector3d dipole;
  std::shared_ptr<DotVector> dots;
  Eigen::Array3d grid_spacing;
  PointSetup()
      : damping(std::make_pair(0.0, 0.0)),
        dipole(Eigen::Vector3d::Zero()),
        dots(std::make_shared<DotVector>(DotVector(
            {QuantumDot(Eigen::Vector3d(-0.5, -0.5, 0), 1.0, damping, dipole),
             QuantumDot(Eigen::Vector3d(-0.2, 1.3, 0), 1.0, damping, dipole),
             QuantumDot(Eigen::Vector3d(0.2, 1.3, 0), 1.0, damping, dipole),
             QuantumDot(Eigen::Vector3d(0.4, -0.7, 0), 1.0, damping, dipole),
             QuantumDot(Eigen::Vector3d(0.6, 0.6, 0), 1.0, damping, dipole),
             QuantumDot(Eigen::Vector3d(0.6, 1.6, 0), 1.0, damping, dipole),
             QuantumDot(Eigen::Vector3d(0.6, 0.1, 0), 1.0, damping, dipole)}))),
        grid_spacing(Eigen::Vector3d(1, 1, 1)){};
};

BOOST_FIXTURE_TEST_CASE(Construction, PointSetup)
{
  AIM::Grid grid(grid_spacing, dots);

  BOOST_CHECK_EQUAL(grid.dimensions(0), 2);
  BOOST_CHECK_EQUAL(grid.dimensions(1), 3);
  BOOST_CHECK_EQUAL(grid.dimensions(2), 1);

  BOOST_CHECK_CLOSE(grid.max_diagonal, std::sqrt(14.0), 1e-16);

  BOOST_CHECK_EQUAL(grid.boxes[0].second - grid.boxes[0].first, 1);
  BOOST_CHECK_EQUAL(grid.boxes[1].second - grid.boxes[1].first, 1);
  BOOST_CHECK_EQUAL(grid.boxes[2].second - grid.boxes[2].first, 0);
  BOOST_CHECK_EQUAL(grid.boxes[3].second - grid.boxes[3].first, 2);
  BOOST_CHECK_EQUAL(grid.boxes[4].second - grid.boxes[4].first, 1);
  BOOST_CHECK_EQUAL(grid.boxes[5].second - grid.boxes[5].first, 2);
}

BOOST_FIXTURE_TEST_CASE(PointSort, PointSetup)
{
  AIM::Grid grid(grid_spacing, dots);

  auto grid_idx =
      [&](const Eigen::Vector3d &p) {
        return grid.coord_to_idx(grid.grid_coordinate(p));
      };

  BOOST_CHECK_EQUAL(grid_idx(dots->at(0).position()), 0);
  BOOST_CHECK_EQUAL(grid_idx(dots->at(1).position()), 1);
  BOOST_CHECK_EQUAL(grid_idx(dots->at(2).position()), 3);
  BOOST_CHECK_EQUAL(grid_idx(dots->at(3).position()), 3);
  BOOST_CHECK_EQUAL(grid_idx(dots->at(4).position()), 4);
  BOOST_CHECK_EQUAL(grid_idx(dots->at(5).position()), 5);
  BOOST_CHECK_EQUAL(grid_idx(dots->at(6).position()), 5);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
