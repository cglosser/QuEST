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

  fftw_execute_dft(aim.spatial_transforms.forward,
                   reinterpret_cast<fftw_complex *>(nums.data()),
                   reinterpret_cast<fftw_complex *>(nums.data()));
  fftw_execute_dft(aim.spatial_transforms.backward,
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

BOOST_AUTO_TEST_SUITE_END()  // StaticPropagation

BOOST_AUTO_TEST_SUITE_END()  // AIM
