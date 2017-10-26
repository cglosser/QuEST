#include <boost/test/unit_test.hpp>
#include <cmath>

#include "interactions/AIM/aim_interaction.h"

struct PARAMETERS {
  double c, dt;
  int interpolation_order, expansion_order;
  PARAMETERS() : c(1), dt(1), interpolation_order(4), expansion_order(0){};
};

BOOST_AUTO_TEST_SUITE(AIM)

BOOST_AUTO_TEST_CASE(GAUSSIAN_PROPAGATION)
{
  const int interp_order = 4;
  Eigen::Vector3i num_boxes(4, 4, 4);
  Grid grid(Eigen::Vector3d(16, 16, 16), num_boxes);
  AIM::AimInteraction aim(interp_order, grid, AIM::normalization::unit);
  auto circulant_shape = grid.circulant_shape(1, 1);

  // Set the source signal
  double mean = grid.max_diagonal / 2.0;
  double sd = mean / 5;
  auto src = [=](const double t) {
    return std::exp(std::pow((t - mean) / sd, 2) / -2.0);
  };

  // Source signal assumed to radiate from 0th grid point
  for(int t = 0; t < circulant_shape[0]; ++t) {
    aim.source_table[t][0][0][0] = src(t);
    fftw_execute_dft(
        aim.spatial_transforms.forward,
        reinterpret_cast<fftw_complex *>(&aim.source_table[t][0][0][0]),
        reinterpret_cast<fftw_complex *>(&aim.source_table[t][0][0][0]));
  }

  // Check against analytic retardation of source
  for(int t = 0; t < circulant_shape[0]; ++t) {
    auto retarded_signal = aim.evaluate(t).real();
    for(auto b = 1u; b < grid.num_boxes; ++b) {
      const double delay =
          (grid.spatial_coord_of_box(b) - grid.spatial_coord_of_box(0)).norm();
      // Must use BOOST_CHECK_SMALL because some retarded_signal values may be
      // zero (as the signal hasn't arrived yet)
      BOOST_CHECK_SMALL(src(t - delay) - retarded_signal(b), 1e-5);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(GAUSSIAN_POINT_PROPAGATION, PARAMETERS)
{
  Eigen::Array3i num_boxes(4, 4, 4);
  Eigen::Array3d spacing(Eigen::Array3d(1, 1, 1) * c * dt);

  // Place one QD *on* the most-separated grid points
  std::shared_ptr<DotVector> dots = std::make_shared<DotVector>(
      DotVector{QuantumDot(Eigen::Vector3d::Zero()),
                QuantumDot(spacing * num_boxes.cast<double>())});

  Grid grid(spacing, dots, expansion_order);
  BOOST_CHECK_EQUAL(
      (dots->at(1).position() - grid.spatial_coord_of_box(grid.num_boxes - 1))
          .norm(),
      0);

  auto expansions =
      LeastSquaresExpansionSolver::get_expansions(expansion_order, grid, *dots);

  const int num_steps = 64;

  auto src = [=](const double t) {
    int i = t / dt;
    double arg = (i - num_steps / 2.0) / (num_steps / 5.0);
    return gaussian(arg);
  };

  // Set up and pre-fill the source particle in a History table
  std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history =
      std::make_shared<Integrator::History<Eigen::Vector2cd>>(dots->size(), 10,
                                                              num_steps);
  history->fill(Eigen::Vector2cd::Zero());
  for(int i = -10; i < num_steps; ++i) {
    history->array[0][i][0](RHO_01) = src(i * dt);
  }

  AIM::AimInteraction aim(dots, history, nullptr, interpolation_order, c, dt,
                          grid, expansions, AIM::normalization::unit);

  for(int i = 0; i < num_steps; ++i) {
    aim.fill_source_table(i);
    aim.evaluate(i);
  }
}

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

BOOST_AUTO_TEST_CASE(VectorFourierTransforms)
{
  Eigen::Vector3i num_boxes(4, 4, 4);
  Grid grid(unit_spacing, num_boxes);
  auto expansions =
      LeastSquaresExpansionSolver::get_expansions(expansion_order, grid, *dots);
  AIM::AimInteraction aim(dots, history, propagator, interp_order, c0, dt, grid,
                          expansions, AIM::normalization::unit);
  auto circulant_shape = grid.circulant_shape(c0, dt);

  std::fill(aim.source_table.data(),
            aim.source_table.data() + aim.source_table.num_elements(),
            cmplx(0, 0));

  for(int t = 0; t < circulant_shape[0]; ++t) {
    int i = 1;
    for(int x = 0; x < circulant_shape[1] / 2; ++x) {
      for(int y = 0; y < circulant_shape[2] / 2; ++y) {
        for(int z = 0; z < circulant_shape[3] / 2; ++z) {
          aim.source_table[t][x][y][z] = i++;
        }
      }
    }
  }

  fftw_complex *ptr = reinterpret_cast<fftw_complex *>(aim.source_table.data());
  fftw_execute_dft(aim.spatial_transforms.forward, ptr, ptr);

  std::vector<cmplx> test_fft = {
#include "range_16_fft.dat"
  };  // This is cheeky. Don't do it often.

  // Computed with 30-digit precision in Mathematica
  std::vector<cmplx> test_int = {
#include "range_16_int.dat"
  };

  // The spatial transform should transform the first "time block"...
  for(auto j = 0u; j < test_fft.size(); ++j) {
    BOOST_CHECK_SMALL(std::abs(test_fft.at(j) - aim.source_table.data()[j]),
                      1e-12);
  }

  // ...and leave the second one alone.
  for(auto j = 0u; j < test_int.size(); ++j) {
    BOOST_CHECK_EQUAL(test_int.at(j), (&aim.source_table[1][0][0][0])[j]);
  }
}

BOOST_AUTO_TEST_SUITE_END()  // AIM_Fourier_transforms

BOOST_AUTO_TEST_SUITE_END()  // AIM
