#include <boost/test/unit_test.hpp>
#include <cmath>

#include "interactions/AIM/aim_interaction.h"

BOOST_AUTO_TEST_SUITE(AIM)

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
  auto expansions =
      LeastSquaresExpansionSolver::get(expansion_order, grid, *dots);

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
  auto expansions =
      LeastSquaresExpansionSolver::get(expansion_order, grid, *dots);

  for(int dot = 0; dot < 2; ++dot) {
    for(int pt = 0; pt < 8; ++pt) {
      BOOST_CHECK_EQUAL(expansions[dot][pt].weight, 1.0 / 8);
    }
  }
}

BOOST_AUTO_TEST_CASE(VectorFourierTransforms)
{
  Eigen::Vector3i num_boxes(4, 4, 4);
  Grid grid(unit_spacing, num_boxes);
  auto expansions =
      LeastSquaresExpansionSolver::get(expansion_order, grid, *dots);
  AIM::AimInteraction aim(dots, history, propagator, interp_order, c0, dt, grid,
                          expansions);
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

BOOST_AUTO_TEST_CASE(GAUSSIAN_PROPAGATION)
{
  const double step = 0.5;
  Eigen::Vector3i num_boxes(4, 4, 4);
  Grid grid(unit_spacing, num_boxes);
  auto expansions =
      LeastSquaresExpansionSolver::get(expansion_order, grid, *dots);
  AIM::AimInteraction aim(dots, history, propagator, interp_order, c0, step,
                          grid, expansions);
  auto circulant_shape = grid.circulant_shape(c0, step);

  // Set the source radiator -- presumed to sit on a grid point
  double mean = 3.5, sd = 1;
  for(int t = 0; t < circulant_shape[0]; ++t) {
    double val = std::exp(std::pow((t * step - mean) / sd, 2) / -2.0);
    aim.source_table[t][0][0][0] = val;
    fftw_execute_dft(
        aim.spatial_transforms.forward,
        reinterpret_cast<fftw_complex *>(&aim.source_table[t][0][0][0]),
        reinterpret_cast<fftw_complex *>(&aim.source_table[t][0][0][0]));
  }

  for(int t = 0; t < circulant_shape[0]; ++t) {
    // std::cout << std::exp(std::pow((t * step - mean) / sd, 2) / -2.0) << "
    // 0";
    std::cout << t << " " << t * step << " "
              << aim.evaluate(t).transpose().real() << std::endl;
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

BOOST_AUTO_TEST_SUITE_END()  // StaticPropagation

BOOST_AUTO_TEST_SUITE_END()  // AIM
