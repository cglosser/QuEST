#include <boost/test/unit_test.hpp>
#include <cmath>

#include "interactions/AIM/grid.h"

BOOST_AUTO_TEST_SUITE(GRID)

BOOST_AUTO_TEST_SUITE(EUCLIDEAN_COORDINATES)

BOOST_AUTO_TEST_CASE(COORDINATE_TRANSFORMATIONS_2D)
{
  Eigen::Vector3d spacing(1, 1, 1);
  Eigen::Array3i shape(6, 5, 1), shift(-2, -2, 0);
  AIM::Grid grid(spacing, shape, shift);

  // Check "bottom left" gridpt maps to index 0
  BOOST_CHECK_EQUAL(grid.idx_to_coord(0), Eigen::Vector3i::Zero());
  BOOST_CHECK_EQUAL(grid.coord_to_idx(Eigen::Vector3i::Zero()), 0);
  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(0),
                    (spacing.array() * shift.cast<double>()).matrix());

  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(12), Eigen::Vector3d::Zero());

  // Check that "upper right" gridpt maps to last valid index
  BOOST_CHECK_EQUAL(
      grid.spatial_coord_of_box(grid.num_gridpoints - 1),
      grid.spatial_coord_of_box(0) +
          (spacing.array() * (shape.cast<double>() - 1)).matrix());

  // Check that coord_to_idx and idx_to_coord are inverses
  for(size_t i = 0; i < grid.num_gridpoints; ++i) {
    BOOST_CHECK_EQUAL(grid.coord_to_idx(grid.idx_to_coord(i)), i);
  }
}

BOOST_AUTO_TEST_CASE(COORDINATE_TRANSFORMATIONS_3D)
{
  Eigen::Vector3d spacing(1, 1, 1);
  Eigen::Array3i shape(6, 5, 2), shift(-2, -3, 1);
  AIM::Grid grid(spacing, shape, shift);

  // Check "bottom left" gridpt maps to index 0
  BOOST_CHECK_EQUAL(grid.idx_to_coord(0), Eigen::Vector3i::Zero());
  BOOST_CHECK_EQUAL(grid.coord_to_idx(Eigen::Vector3i::Zero()), 0);
  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(0),
                    (spacing.array() * shift.cast<double>()).matrix());

  // Check that "upper right" gridpt maps to last valid index
  BOOST_CHECK_EQUAL(
      grid.spatial_coord_of_box(grid.num_gridpoints - 1),
      grid.spatial_coord_of_box(0) +
          (spacing.array() * (shape.cast<double>() - 1)).matrix());

  // Check that coord_to_idx and idx_to_coord are inverses
  for(size_t i = 0; i < grid.num_gridpoints; ++i) {
    BOOST_CHECK_EQUAL(grid.coord_to_idx(grid.idx_to_coord(i)), i);
  }
}

BOOST_AUTO_TEST_CASE(NONINTEGER_SHIFT)
{
  Eigen::Vector3d spacing(0.2, 0.2, 0.2);
  Eigen::Array3i shape(10, 10, 10), shift(-5, -5, -5);
  AIM::Grid grid(spacing, shape, shift);

  // Origin is invariant
  BOOST_CHECK_EQUAL(
      grid.grid_coordinate(Eigen::Vector3d::Zero()).cast<double>(),
      Eigen::Vector3d::Zero());

  // First orthant
  BOOST_CHECK_EQUAL(
      grid.grid_coordinate(Eigen::Vector3d(1, 2, 2.1)).cast<double>(),
      Eigen::Vector3d(5, 10, 10));

  // Opposite orthant
  BOOST_CHECK_EQUAL(
      grid.grid_coordinate(Eigen::Vector3d(-2.3, -3.7, -1.4)).cast<double>(),
      Eigen::Vector3d(-11, -18, -6));

  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(0),
                    (shift.cast<double>() * spacing.array()).matrix());
}

BOOST_AUTO_TEST_SUITE_END()  // EUCLIDEAN_COORDINATES

BOOST_AUTO_TEST_CASE(CIRCULANT_MATRIX_SHAPE)
{
  Eigen::Vector3d spacing(1, 1, 1);
  Eigen::Vector3i shape(20, 20, 20);
  AIM::Grid grid(spacing, shape);

  auto dims = grid.circulant_shape(1, 1);
  double diag = shape.cast<double>().norm();
  BOOST_CHECK_EQUAL(dims[0], std::ceil(diag));
  BOOST_CHECK_EQUAL(dims[1], 2 * shape[0]);
  BOOST_CHECK_EQUAL(dims[2], 2 * shape[1]);
  BOOST_CHECK_EQUAL(dims[3], 2 * shape[2]);
}

struct PARAMETERS {
  Eigen::Array3d spacing;
  double c, dt;
  int expansion_order;

  std::shared_ptr<DotVector> dots;

  PARAMETERS()
      : spacing(1, 1, 1),
        c(1),
        dt(1),
        expansion_order(1),
        dots(std::make_shared<DotVector>()){};
};

BOOST_FIXTURE_TEST_SUITE(DOTS, PARAMETERS)

BOOST_AUTO_TEST_CASE(EXPANSION_COORDINATES)
{
  dots->push_back(QuantumDot(Eigen::Vector3d(0, -3, -3)));
  dots->push_back(QuantumDot(Eigen::Vector3d(0, 3, 3)));
  AIM::Grid grid(spacing, dots, expansion_order);

  BOOST_CHECK_EQUAL(grid.coord_to_idx(grid.idx_to_coord(0)), 0);

  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(0), Eigen::Vector3d(0, -3, -3));
  BOOST_CHECK_EQUAL(
      grid.spatial_coord_of_box(grid.num_gridpoints - 1),
      Eigen::Vector3d(1, 4, 4));  // <1, 4, 4> due to the expansion order

  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(
                        grid.associated_grid_index(dots->at(0).position())),
                    dots->at(0).position());
  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(
                        grid.associated_grid_index(dots->at(1).position())),
                    dots->at(1).position());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_CASE(EXPANSION_DISTANCE)
{
  AIM::Grid grid1(Eigen::Vector3d(1, 1, 1), Eigen::Array3i(10, 10, 10), Eigen::Vector3i::Zero(), 1);

  BOOST_CHECK_EQUAL(grid1.expansion_distance(0, 0), 0);
  BOOST_CHECK_EQUAL(grid1.expansion_distance(0, 1), 0);
  BOOST_CHECK_EQUAL(grid1.expansion_distance(0, 2), 1);
  BOOST_CHECK_EQUAL(grid1.expansion_distance(0, 8), 7);
  BOOST_CHECK_EQUAL(grid1.expansion_distance(0, 10), 0);
  BOOST_CHECK_EQUAL(grid1.expansion_distance(0, 80), 7);
  BOOST_CHECK_EQUAL(grid1.expansion_distance(0, 88), 7);
  BOOST_CHECK_EQUAL(grid1.expansion_distance(0, 888), 7);

  AIM::Grid grid2(Eigen::Vector3d(1, 1, 1), Eigen::Array3i(5, 5, 5), Eigen::Vector3i(-1, -1, -1), 2);

  BOOST_CHECK_EQUAL(grid2.expansion_distance(31, 31), 0);
  BOOST_CHECK_EQUAL(grid2.expansion_distance(31, 32), 0);
  BOOST_CHECK_EQUAL(grid2.expansion_distance(31, 33), 0);
  BOOST_CHECK_EQUAL(grid2.expansion_distance(31, 99), 1);
}

BOOST_AUTO_TEST_CASE(NEARFIELD_PAIRS)
{
  auto dots = std::make_shared<DotVector>();
  dots->push_back(QuantumDot(Eigen::Vector3d(0.5, 0.5, 0.5), 0, {0, 0}, Eigen::Vector3d(0,0,1)));
  dots->push_back(QuantumDot(Eigen::Vector3d(0.5, 0.5, 1.5), 0, {0, 0}, Eigen::Vector3d(0,0,1)));
  dots->push_back(QuantumDot(Eigen::Vector3d(0.5, 0.5, 7.5), 0, {0, 0}, Eigen::Vector3d(0,0,1)));

  AIM::Grid grid(Eigen::Vector3d(1, 1, 1), dots, 1);
  const auto nf = grid.nearfield_pairs(2);

  BOOST_CHECK_EQUAL(nf.at(0).first, 0);
  BOOST_CHECK_EQUAL(nf.at(0).second, 0);

  BOOST_CHECK_EQUAL(nf.at(1).first, 0);
  BOOST_CHECK_EQUAL(nf.at(1).second, 1);

  BOOST_CHECK_EQUAL(nf.at(2).first, 1);
  BOOST_CHECK_EQUAL(nf.at(2).second, 1);

  BOOST_CHECK_EQUAL(nf.at(3).first, 7);
  BOOST_CHECK_EQUAL(nf.at(3).second, 7);

  BOOST_CHECK_EQUAL(nf.size(), 4);

}

BOOST_AUTO_TEST_SUITE_END()  // GRID
