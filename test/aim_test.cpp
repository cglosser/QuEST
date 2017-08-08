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

  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(0), Eigen::Vector3d(-1, -1, 0));
  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(1), Eigen::Vector3d(0, -1, 0));
  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(2), Eigen::Vector3d(-1, 0, 0));
  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(3), Eigen::Vector3d(0, 0, 0));
  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(4), Eigen::Vector3d(-1, 1, 0));
  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(5), Eigen::Vector3d(0, 1, 0));
}

BOOST_FIXTURE_TEST_CASE(PointSort, PointSetup)
{
  AIM::Grid grid(grid_spacing, dots);

  auto grid_idx = [&](const Eigen::Vector3d &p) {
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

struct Universe {
  double c, dt;
  std::shared_ptr<DotVector> dots;
  Eigen::Array3d grid_spacing;
  Universe()
      : c(1),
        dt(1),
        dots(std::make_shared<DotVector>()),
        grid_spacing(1.2, 1.2, 1.2)
  {
    for(int x = 0; x < 6; ++x) {
      for(int y = 0; y < 6; ++y) {
        dots->push_back(QuantumDot(1.2 * Eigen::Vector3d(x, y, 0), 0,
                                   std::make_pair(10.0, 20.0),
                                   Eigen::Vector3d(0, 0, 0)));
      }
    }
  };
};

BOOST_FIXTURE_TEST_CASE(two_d_matrix_elements, Universe)
{
  AIM::Grid grid(grid_spacing, dots);
  AIM::AimInteraction aim(dots, grid_spacing, 3, 1, 1);

  auto basis_vals = aim.g_matrix_row(1);

  const double v1 = 84./125, v2 = 2 * (520 - 361 * sqrt(2.0)) / 125;

  for(size_t i = 0; i < basis_vals.size(); ++i) {
    switch(i) {
      case 1: BOOST_CHECK_CLOSE(basis_vals[i], v1, 1e-12); break;
      case 6: BOOST_CHECK_CLOSE(basis_vals[i], v1, 1e-12); break;
      case 7: BOOST_CHECK_CLOSE(basis_vals[i], v2, 1e-12); break;
      case 36: BOOST_CHECK_CLOSE(basis_vals[i], v1, 1e-12); break;
      case 37: BOOST_CHECK_CLOSE(basis_vals[i], v2, 1e-12); break;
      case 42: BOOST_CHECK_CLOSE(basis_vals[i], v2, 1e-12); break;
      default: BOOST_CHECK_EQUAL(basis_vals[i], 0); break;
    }
  }
}

BOOST_FIXTURE_TEST_CASE(diagonal_representation, Universe)
{
  AIM::AimInteraction aim(dots, grid_spacing, 3, 1, 1);

  std::cout << aim.fourier_table.shape()[0] << " "
            << aim.fourier_table.shape()[1] << std::endl;

  for(int i = 0; i < 24; ++i) {
    std::cout << aim.fourier_table[1][i] << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
