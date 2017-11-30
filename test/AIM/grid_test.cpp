#include <boost/test/unit_test.hpp>
#include <cmath>

#include "interactions/AIM/grid.h"

BOOST_AUTO_TEST_SUITE(GridTest)

BOOST_AUTO_TEST_CASE(DummyConstructor)
{
  AIM::Grid grid(Eigen::Vector3d(1, 1, 1), Eigen::Vector3i(30, 30, 30));
  BOOST_CHECK_EQUAL(grid.num_gridpoints, std::pow(30, 3));
}

BOOST_AUTO_TEST_CASE(CIRCULANT_MATRIX_SHAPE)
{
  AIM::Grid grid(Eigen::Vector3d(1, 1, 1), Eigen::Vector3i(20, 20, 20));

  auto dims = grid.circulant_shape(1, 1);
  BOOST_CHECK_EQUAL(dims[0], 35);
  BOOST_CHECK_EQUAL(dims[1], 40);
  BOOST_CHECK_EQUAL(dims[2], 40);
  BOOST_CHECK_EQUAL(dims[3], 40);
}

BOOST_AUTO_TEST_SUITE(SingleDotExpansions)

struct OffsetDot {
  std::shared_ptr<DotVector> dot;
  Eigen::Vector3d unit_spacing;
  int max_order;
  OffsetDot()
      : dot(std::make_shared<DotVector>()),
        unit_spacing(Eigen::Vector3d(1, 1, 1)),
        max_order(50)
  {
    dot->push_back(QuantumDot(Eigen::Vector3d(0.5, 0.5, 0.5), 0, {0.0, 0.0},
                              Eigen::Vector3d(0, 0, 0)));
  }
};

BOOST_FIXTURE_TEST_CASE(GRID_SUFFICIENT_FOR_EXPANSIONS, OffsetDot)
{
  // Check that the grid contains enough points for the expansion of the qdot

  for(int order = 0; order < max_order; ++order) {
    AIM::Grid grid(unit_spacing, dot, order);
    BOOST_CHECK_EQUAL(grid.num_gridpoints, std::pow(order + 1, 3));
  }
}

BOOST_FIXTURE_TEST_CASE(AssocBoxPoint, OffsetDot)
{
  // Check that the QD always corresponds to the gridpoint located at (0,0,0)
  // in space

  for(int order = 0; order < max_order; ++order) {
    AIM::Grid grid(unit_spacing, dot, order);
    auto idx = grid.coord_to_idx(grid.grid_coordinate(dot->front().position()));
    auto spatial = grid.spatial_coord_of_box(idx);

    BOOST_CHECK_EQUAL(spatial(0), 0);
    BOOST_CHECK_EQUAL(spatial(1), 0);
    BOOST_CHECK_EQUAL(spatial(2), 0);
  }
}

BOOST_FIXTURE_TEST_CASE(In_middle_box_for_even_expansions, OffsetDot)
{
  // Check that the quantum dot always lies in the "middle" box when expanded
  // to an even order

  for(int order = 0; order < max_order; order += 2) {
    AIM::Grid grid(unit_spacing, dot, order);
    auto idx = grid.coord_to_idx(grid.grid_coordinate(dot->front().position()));
    BOOST_CHECK_EQUAL(idx, grid.num_gridpoints / 2);
  }
}

BOOST_FIXTURE_TEST_CASE(ExpansionIndices, OffsetDot)
{
  // Check that a single dot expands "into" [0, num_gridpoints) available gridpoints

  for(int order = 0; order < max_order; ++order) {
    AIM::Grid grid(unit_spacing, dot, order);
    auto expansion_indices =
        grid.expansion_box_indices(dot->front().position(), order);

    BOOST_CHECK_EQUAL(expansion_indices.size(), std::pow(order + 1, 3));

    // This isn't an ideal check as it depends on the convention in
    // expansion_box_indices(...) (so the indices might have a weird ordering),
    // but it at least provides some bounds checking and can warn if things go
    // horribly wrong
    std::sort(expansion_indices.begin(), expansion_indices.end());
    for(size_t i = 0; i < expansion_indices.size(); ++i) {
      BOOST_CHECK_EQUAL(expansion_indices.at(i), i);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()  // SingleDotExpansions

BOOST_AUTO_TEST_SUITE(TwoDotExpansions)

struct TwoOffsetDots {
  std::shared_ptr<DotVector> dots;
  int max_order;
  TwoOffsetDots() : dots(std::make_shared<DotVector>()), max_order(50)
  {
    dots->push_back(QuantumDot(Eigen::Vector3d(0.5, 0.5, 0.5), 0, {0.0, 0.0},
                               Eigen::Vector3d(0, 0, 0)));
    dots->push_back(QuantumDot(Eigen::Vector3d(0.5, 0.5, 1.5), 0, {0.0, 0.0},
                               Eigen::Vector3d(0, 0, 0)));
  }
};


BOOST_FIXTURE_TEST_CASE(NumberOfGridPts, TwoOffsetDots)
{
  for(int order = 0; order < max_order; ++order) {
    AIM::Grid grid(Eigen::Vector3d(1, 1, 1), dots, order);
    BOOST_CHECK_EQUAL(grid.num_gridpoints, std::pow(order + 1, 2) * (order + 2));
  }
}

BOOST_AUTO_TEST_SUITE_END()  // TwoDotExpansions

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
  AIM::Grid grid(grid_spacing, dots, 0);

  auto boxes = grid.box_contents_map(dots);

  BOOST_CHECK_EQUAL(grid.dimensions(0), 2);
  BOOST_CHECK_EQUAL(grid.dimensions(1), 3);
  BOOST_CHECK_EQUAL(grid.dimensions(2), 1);

  BOOST_CHECK_CLOSE(grid.max_diagonal, std::sqrt(14.0), 1e-16);

  BOOST_CHECK_EQUAL(boxes[0].second - boxes[0].first, 1);
  BOOST_CHECK_EQUAL(boxes[1].second - boxes[1].first, 1);
  BOOST_CHECK_EQUAL(boxes[2].second - boxes[2].first, 0);
  BOOST_CHECK_EQUAL(boxes[3].second - boxes[3].first, 2);
  BOOST_CHECK_EQUAL(boxes[4].second - boxes[4].first, 1);
  BOOST_CHECK_EQUAL(boxes[5].second - boxes[5].first, 2);

  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(0), Eigen::Vector3d(-1, -1, 0));
  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(1), Eigen::Vector3d(0, -1, 0));
  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(2), Eigen::Vector3d(-1, 0, 0));
  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(3), Eigen::Vector3d(0, 0, 0));
  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(4), Eigen::Vector3d(-1, 1, 0));
  BOOST_CHECK_EQUAL(grid.spatial_coord_of_box(5), Eigen::Vector3d(0, 1, 0));
}

BOOST_FIXTURE_TEST_CASE(PointSort, PointSetup)
{
  AIM::Grid grid(grid_spacing, dots, 0);

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

BOOST_AUTO_TEST_SUITE_END()  // TwoDimensions

BOOST_AUTO_TEST_SUITE_END()  // GridTest
