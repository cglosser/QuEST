#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iomanip>

#include "interactions/AIM/aim_interaction.h"

struct PARAMETERS {
  double c, dt;
  int interpolation_order, expansion_order;
  PARAMETERS() : c(1), dt(1), interpolation_order(4), expansion_order(0){};
};

BOOST_AUTO_TEST_SUITE(AIM)

BOOST_FIXTURE_TEST_CASE(GAUSSIAN_POINT_PROPAGATION, PARAMETERS)
{
  Eigen::Array3i num_boxes(4, 4, 4);
  Eigen::Array3d spacing(Eigen::Array3d(1, 1, 1) * c * dt);

  // Place one QD *on* the most-separated grid points
  std::shared_ptr<DotVector> dots = std::make_shared<DotVector>(
      DotVector{QuantumDot(Eigen::Vector3d::Zero(), Eigen::Vector3d(0, 0, 1)),
                QuantumDot(spacing * num_boxes.cast<double>(),
                           Eigen::Vector3d(0, 0, 1))});

  Grid grid(spacing, dots, expansion_order);
  BOOST_CHECK_EQUAL(
      (dots->at(1).position() - grid.spatial_coord_of_box(grid.num_boxes - 1))
          .norm(),
      0);

  auto expansions = Expansions::LeastSquaresExpansionSolver::get_expansions(
      expansion_order, grid, *dots);

  const int num_steps = 256;

  auto src = [=](const double t) {
    const double total_time = num_steps * dt;
    double arg = (t - total_time / 2.0) / (total_time / 6.0);
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
                          grid, expansions, AIM::Expansions::Identity,
                          AIM::normalization::unit);

  const double delay =
      (dots->at(1).position() - dots->at(0).position()).norm() / c;

  double max_error = 0;
  for(int i = 0; i < num_steps; ++i) {
    auto x = aim.evaluate(i);

    if(i > aim.max_transit_steps) {
      // Wait until observer is fully (i.e. has a complete interpolation) inside
      // of light cone
      BOOST_CHECK_CLOSE(x(1).real(), src(i * dt - delay), 1e-7);

      auto relative_error =
          std::abs((src(i * dt - delay) - x(1).real()) / src(i * dt - delay));
      max_error = std::max(max_error, relative_error);
    }
  }
  BOOST_TEST_MESSAGE("Maximum relative on-grid error: " << max_error);
}

BOOST_FIXTURE_TEST_CASE(GRAD_DIV, PARAMETERS)
{
  Eigen::Array3i num_boxes(4, 4, 4);
  Eigen::Array3d spacing(Eigen::Array3d(1, 1, 1) * c * dt);

  // Place one QD *on* the most-separated grid points
  const double h = 0.5;
  std::shared_ptr<DotVector> dots = std::make_shared<DotVector>(DotVector{
      QuantumDot(Eigen::Vector3d(h, h, h), Eigen::Vector3d(0, 0, 1)),
      QuantumDot(
          spacing * (num_boxes.cast<double>() + Eigen::Array3d(h, h, h)),
          Eigen::Vector3d(0, 0, 1))});

  Grid grid(spacing, dots, 1);
  auto expansions = Expansions::LeastSquaresExpansionSolver::get_expansions(
      1, grid, *dots);

  const int num_steps = 256;

  auto src = [=](const double t) {
    const double total_time = num_steps * dt;
    double arg = (t - total_time / 2.0) / (total_time / 6.0);
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
                          grid, expansions, AIM::Expansions::GradDiv,
                          AIM::normalization::unit);

  const double delay =
      (dots->at(1).position() - dots->at(0).position()).norm() / c;

  std::cout << std::scientific;
  std::cout << "START HERE" << std::endl;
  for(int i = 0; i < num_steps; ++i) {
    auto x = aim.evaluate(i);
    //std::cout << x.transpose().real() << std::endl;
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
  auto expansions = Expansions::LeastSquaresExpansionSolver::get_expansions(
      expansion_order, grid, *dots);
  AIM::AimInteraction aim(dots, history, propagator, interp_order, c0, dt, grid,
                          expansions, AIM::Expansions::Identity,
                          AIM::normalization::unit);
  auto circulant_shape = grid.circulant_shape(c0, dt);

  std::fill(aim.source_table.data(),
            aim.source_table.data() + aim.source_table.num_elements(),
            cmplx(0, 0));

  for(int t = 0; t < circulant_shape[0]; ++t) {
    int i = 1;
    for(int x = 0; x < circulant_shape[1] / 2; ++x) {
      for(int y = 0; y < circulant_shape[2] / 2; ++y) {
        for(int z = 0; z < circulant_shape[3] / 2; ++z) {
          Eigen::Map<Eigen::Vector3cd> m(&aim.source_table[t][x][y][z][0]);
          m = Eigen::Vector3cd(i, -2 * i, 0);
          i++;
        }
      }
    }
  }

  fftw_complex *ptr = reinterpret_cast<fftw_complex *>(aim.source_table.data());
  fftw_execute_dft(aim.spatial_vector_transforms.forward, ptr, ptr);

  std::vector<cmplx> test_fft = {
#include "range_16_fft.dat"
  };  // Computed with 30-digit precision in Mathematica

  std::vector<cmplx> test_int = {
#include "range_16_int.dat"
  };

  // The spatial transform should transform the first "time block"...
  Eigen::Map<Eigen::Array3Xcd> vecs(aim.source_table.data(), 3,
                                    (2 * num_boxes).prod());
  for(auto j = 0u; j < test_fft.size(); ++j) {
    BOOST_CHECK_SMALL(std::norm(test_fft.at(j) - vecs(0, j)), 1e-16);
    BOOST_CHECK_SMALL(std::norm(-2.0 * test_fft.at(j) - vecs(1, j)), 1e-16);
    BOOST_CHECK_SMALL(std::norm(vecs(2, j)), 1e-16);
  }

  // ...and leave the second one alone.
  new(&vecs) Eigen::Map<Eigen::Array3Xcd>(&aim.source_table[1][0][0][0][0], 3,
                                          (2 * num_boxes).prod());
  for(auto j = 0u; j < test_int.size(); ++j) {
    BOOST_CHECK_EQUAL(test_int.at(j), vecs(0, j));
    BOOST_CHECK_EQUAL(-2.0 * test_int.at(j), vecs(1, j));
    BOOST_CHECK_EQUAL(0.0, vecs(2, j));
  }
}

BOOST_AUTO_TEST_SUITE_END()  // AIM_Fourier_transforms

BOOST_AUTO_TEST_SUITE_END()  // AIM
