#include <boost/test/unit_test.hpp>
#include <cmath>

#include "../src/interactions/aim_interaction.h"

BOOST_AUTO_TEST_SUITE(AIM)

BOOST_AUTO_TEST_SUITE(GridTest)

BOOST_AUTO_TEST_CASE(DummyConstructor)
{
  Grid grid(Eigen::Vector3d(1, 1, 1), Eigen::Vector3i(30, 30, 30));
  BOOST_CHECK_EQUAL(grid.num_boxes, std::pow(30, 3));
}

BOOST_AUTO_TEST_CASE(FourierIndices)
{
  Grid grid(Eigen::Vector3d(1, 1, 1), Eigen::Vector3i(3, 4, 5));

  BOOST_CHECK_EQUAL(grid.toeplitz_matrix_size(), 12);

  for(int i = 0; i < 12; ++i) {
    BOOST_CHECK_EQUAL(grid.fourier_idx(0, i), i);
  }

  for(int i = 0; i < 12; ++i) {
    for(int j = 0; j < 12; ++j) {
      BOOST_CHECK_EQUAL(grid.fourier_idx(i, j), grid.fourier_idx(j, i));
    }
  }

  BOOST_CHECK_EQUAL(grid.fourier_idx(1, 0), 1);
  BOOST_CHECK_EQUAL(grid.fourier_idx(1, 1), 0);
  BOOST_CHECK_EQUAL(grid.fourier_idx(1, 2), 1);
  BOOST_CHECK_EQUAL(grid.fourier_idx(1, 3), 2);
  BOOST_CHECK_EQUAL(grid.fourier_idx(1, 4), 5);
  BOOST_CHECK_EQUAL(grid.fourier_idx(1, 5), 4);
  BOOST_CHECK_EQUAL(grid.fourier_idx(1, 6), 5);
  BOOST_CHECK_EQUAL(grid.fourier_idx(1, 7), 6);
  BOOST_CHECK_EQUAL(grid.fourier_idx(1, 8), 9);
  BOOST_CHECK_EQUAL(grid.fourier_idx(1, 9), 8);
  BOOST_CHECK_EQUAL(grid.fourier_idx(1, 10), 9);
  BOOST_CHECK_EQUAL(grid.fourier_idx(1, 11), 10);
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

BOOST_FIXTURE_TEST_CASE(NumberOfGridPts, OffsetDot)
{
  // Check that the grid contains enough points for the expansion of the qdot

  for(int order = 0; order < max_order; ++order) {
    Grid grid(unit_spacing, dot, order);
    BOOST_CHECK_EQUAL(grid.num_boxes, std::pow(order + 1, 3));
  }
}

BOOST_FIXTURE_TEST_CASE(AssocBoxPoint, OffsetDot)
{
  // Check that the QD always corresponds to the gridpoint located at (0,0,0)
  // in space

  for(int order = 0; order < max_order; ++order) {
    Grid grid(unit_spacing, dot, order);
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
    Grid grid(unit_spacing, dot, order);
    auto idx = grid.coord_to_idx(grid.grid_coordinate(dot->front().position()));
    BOOST_CHECK_EQUAL(idx, grid.num_boxes / 2);
  }
}

BOOST_FIXTURE_TEST_CASE(ExpansionIndices, OffsetDot)
{
  // Check that a single dot expands "into" [0, num_boxes) available gridpoints

  for(int order = 0; order < max_order; ++order) {
    Grid grid(unit_spacing, dot, order);
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

BOOST_AUTO_TEST_SUITE_END()

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
    Grid grid(Eigen::Vector3d(1, 1, 1), dots, order);
    BOOST_CHECK_EQUAL(grid.num_boxes, std::pow(order + 1, 2) * (order + 2));
  }
}

BOOST_AUTO_TEST_SUITE_END()

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

struct DummyPropagation {
  std::shared_ptr<DotVector> dots;
  std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history;
  std::shared_ptr<Propagation::RotatingFramePropagator> propagator;
  int interp_order, expansion_order;
  double c0, dt;
  Eigen::Vector3d unit_spacing;
  DummyPropagation()
      : dots(std::make_shared<DotVector>()),
        history(nullptr),
        propagator(nullptr),
        interp_order(3),
        expansion_order(1),
        c0(1),
        dt(1),
        unit_spacing(1, 1, 1){};
};

BOOST_FIXTURE_TEST_SUITE(AIM_Fourier_transforms, DummyPropagation)

// Checks a couple properties of the arrays used to hold Fourier transform
// data. Does not account for *anything* relating to propagation, therefore
// those data members have been initialized to their appropriate null value.

BOOST_AUTO_TEST_CASE(OnePointExpansion)
{
  dots->push_back(QuantumDot(Eigen::Vector3d(0.5, 0.5, 0.5), 0, {0.0, 0.0},
                             Eigen::Vector3d(0, 0, 0)));
  Grid grid(unit_spacing, dots, expansion_order);
  AIM::AimInteraction aim(dots, history, propagator, interp_order, c0, dt, grid,
                          expansion_order);

  auto expansions = aim.expansions();

  for(int i = 0; i < 8; ++i) {
    BOOST_CHECK_EQUAL(expansions[0][i].weight, 1.0 / 8);
  }
}

BOOST_AUTO_TEST_CASE(TwoPointExpansions)
{
  dots->push_back(QuantumDot(Eigen::Vector3d(0.5, 0.5, 0.5), 0, {0.0, 0.0},
                             Eigen::Vector3d(0, 0, 0)));
  dots->push_back(QuantumDot(Eigen::Vector3d(0.5, 0.5, 10.5), 0, {0.0, 0.0},
                             Eigen::Vector3d(0, 0, 0)));
  Grid grid(unit_spacing, dots, expansion_order);
  AIM::AimInteraction aim(dots, history, propagator, interp_order, c0, dt, grid,
                          expansion_order);

  auto expansions = aim.expansions();

  for(int dot = 0; dot < 2; ++dot) {
    for(int pt = 0; pt < 8; ++pt) {
      BOOST_CHECK_EQUAL(expansions[dot][pt].weight, 1.0 / 8);
    }
  }
}

BOOST_AUTO_TEST_CASE(VectorFourierTransforms)
{
  Eigen::Vector3i num_boxes(20, 20, 20);
  Grid grid(unit_spacing, num_boxes);
  AIM::AimInteraction aim(dots, history, propagator, interp_order, c0, dt, grid,
                          expansion_order);

  std::vector<cmplx> nums(20 * 20 * 40);
  for(size_t i = 0; i < nums.size(); ++i) {
    nums[i] = i;
  }

  fftw_execute_dft(aim.vector_forward_plan,
                   reinterpret_cast<fftw_complex *>(nums.data()),
                   reinterpret_cast<fftw_complex *>(nums.data()));
  fftw_execute_dft(aim.vector_backward_plan,
                   reinterpret_cast<fftw_complex *>(nums.data()),
                   reinterpret_cast<fftw_complex *>(nums.data()));

  for(size_t i = 0; i < nums.size(); ++i) {
    BOOST_CHECK_CLOSE(nums[i].real(), 40.0 * i, 1e-10);
    BOOST_CHECK_SMALL(nums[i].imag(), 1e-10);
  }
}

BOOST_AUTO_TEST_SUITE_END()  // AIM_Fourier_transforms

struct TwoStaticDots {
  double c0, dt;
  int interp_order, expansion_order;
  std::shared_ptr<DotVector> dots;
  std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history;
  std::shared_ptr<Propagation::RotatingFramePropagator> propagator;
  Eigen::Vector3d unit_spacing;
  Grid grid;
  TwoStaticDots()
      : c0(1),
        dt(1),
        interp_order(3),
        expansion_order(1),
        dots(std::make_shared<DotVector>()),
        history(std::make_shared<Integrator::History<Eigen::Vector2cd>>(
            2, 22, 100)),
        propagator(std::make_shared<Propagation::RotatingFramePropagator>(
            1, c0, 1, 1)),
        unit_spacing(Eigen::Vector3d(1, 1, 1))
  {
    dots->push_back(QuantumDot(Eigen::Vector3d(0.5, 0.5, 0.5), 0, {0.0, 0.0},
                               Eigen::Vector3d(0, 0, 0)));
    dots->push_back(QuantumDot(Eigen::Vector3d(0.5, 0.5, 9.5), 0, {0.0, 0.0},
                               Eigen::Vector3d(0, 0, 0)));
    grid = Grid(unit_spacing, dots, expansion_order);
  };
};

BOOST_FIXTURE_TEST_SUITE(StaticPropagation, TwoStaticDots)

BOOST_AUTO_TEST_CASE(Fast_multiply)
{
  std::cout.precision(17);
  std::cout << std::scientific << std::endl;

  AimInteraction aim(dots, history, propagator, interp_order, c0, dt, grid,
                     expansion_order);

  const int zlen = 2 * grid.dimensions(2);

  int i = 0;
  for(int x = 0; x < grid.dimensions(0); ++x) {
    for(int y = 0; y < grid.dimensions(1); ++y) {
      for(int z = 0; z < zlen; ++z) {
        aim.source_table[0][x][y][z] = (z < grid.dimensions(2)) ? i++ : 0;
      }
    }
  }

  cmplx *const src = &aim.source_table[0][0][0][0];
  fftw_execute_dft(aim.vector_forward_plan,
                   reinterpret_cast<fftw_complex *const>(src),
                   reinterpret_cast<fftw_complex *const>(src));

  std::cout << "Hello!" << std::endl;
  auto matmul = aim.fast_multiply(1, 0);

  std::cout << matmul << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()  // StaticPropagation

BOOST_AUTO_TEST_SUITE_END()  // AIM
