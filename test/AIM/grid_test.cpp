#include <boost/test/unit_test.hpp>
#include <cmath>

#include "interactions/AIM/grid.h"

BOOST_AUTO_TEST_SUITE(GRID)

BOOST_AUTO_TEST_CASE(INDEXING_WITH_NEGATIVE_SHIFT)
{
  Eigen::Vector3d spacing(1, 1, 1);
  Eigen::Array3i shape(20, 20, 20), shift(-4, -4, -4);
  AIM::Grid grid(spacing, shape, shift);

  BOOST_CHECK_EQUAL(grid.num_gridpoints, shape.prod());
  BOOST_CHECK_EQUAL(grid.idx_to_coord(0), Eigen::Vector3i::Zero());

  // Check that "bottom left" gridpoint maps to index 0
  int origin = 0;
  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(origin),
                    (shift.cast<double>() * spacing.array()).matrix());

  // Check that "top right" gridpoint maps to index num_pts - 1
  int far_corner = shape.prod() - 1;
  BOOST_CHECK_EQUAL(
      grid.spatial_coord_of_box(far_corner),
      ((shape + shift - 1).cast<double>() * spacing.array()).matrix());

  // Check that coord_to_idx and idx_to_coord are inverses
  for(int i = 0; i < far_corner; i++) {
    BOOST_CHECK_EQUAL(grid.coord_to_idx(grid.idx_to_coord(i)), i);
  }
}

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
  Eigen::Vector3d spacing;
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


BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()  // GRID
