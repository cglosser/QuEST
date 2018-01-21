#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iomanip>

#include "interactions/AIM/aim_interaction.h"

BOOST_AUTO_TEST_SUITE(AIM)

struct PARAMETERS {
  static int interpolation_order, num_steps, num_dots;
  static double c, dt, total_time;

  static Eigen::Array3i num_boxes;
  static Eigen::Array3d spacing;

  int expansion_order;  // Different orders for different test geometries

  std::shared_ptr<DotVector> dots;
  std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history;

  Grid grid;
  Expansions::ExpansionTable expansions;

  PARAMETERS(const int expansion_order, std::shared_ptr<DotVector> dots)
      : expansion_order(expansion_order),
        dots(dots),
        history(std::make_shared<Integrator::History<Eigen::Vector2cd>>(
            num_dots, 10, num_steps)),
        grid(spacing, dots, expansion_order),
        expansions(Expansions::LeastSquaresExpansionSolver::get_expansions(
            expansion_order, grid, *dots)){};

  virtual ~PARAMETERS() = 0;

  double src(const double t) const
  {
    double arg = (t - total_time / 2.0) / (total_time / 6.0);
    return gaussian(arg);
  }

  double dt_src(const double t) const
  {
    const double mu = total_time / 2.0, sigma = total_time / 6.0;
    const double arg = (t - mu) / sigma;
    return -arg * gaussian(arg) / sigma;
  }

  double efie_src(const Eigen::Vector3d &dr, const double t)
  {
    // This is horrible; don't even try to pull it apart yourself.
    // I just did it with Mathematica's symbolics and CForm.
    const double csq = std::pow(c, 2);
    const double r = dr.norm();
    const double x = dr(0), y = dr(1), z = dr(2);
    return (18 * (csq * std::pow(total_time, 2) *
                      (2 * r + c * (total_time - 2 * t)) *
                      (std::pow(x, 2) + std::pow(y, 2)) +
                  8 * r * (3 * r + c * (total_time - 3 * t)) *
                      (2 * c * total_time + 3 * r - 3 * c * t) * (r - z) *
                      (r + z))) /
           (csq * std::exp((9 * std::pow(2 * r + c * (total_time - 2 * t), 2)) /
                           (2 * csq * std::pow(total_time, 2))) *
            std::pow(total_time, 4) * std::pow(r, 3));
  }

  double test_propagation(AIM::AimInteraction &aim,
                          const std::function<double(const double)> &sol,
                          const double toler)
  {
    double max_relative_error = 0;
    for(int i = 0; i < num_steps; ++i) {
      auto x = aim.evaluate(i);

      if(i > grid.max_transit_steps(c, dt) + 2 * interpolation_order) {
        BOOST_CHECK_CLOSE(x(1).real(), sol(i * dt), toler);

        double diff = sol(i * dt) - x(1).real();
        auto relative_error = std::abs(diff) / sol(i * dt);
        max_relative_error = std::max(max_relative_error, relative_error);
      }
    }

    return max_relative_error;
  }
};

// These are common to a large suie of tests and some test parameters (like dot
// positions) might depend on them. To resolve this with a minimum of code
// duplication, these variables have been made STATIC so that they're available
// to all subclasses of PARAMETERS when they're constructed.
PARAMETERS::~PARAMETERS() {}
int PARAMETERS::interpolation_order = 3;
int PARAMETERS::num_steps = 1024;
int PARAMETERS::num_dots = 2;

double PARAMETERS::c = 1;
double PARAMETERS::dt = 1;
double PARAMETERS::total_time = dt * num_steps;

Eigen::Array3i PARAMETERS::num_boxes(8, 8, 8);
Eigen::Array3d PARAMETERS::spacing(Eigen::Array3d(1, 1, 1) * c * dt);

BOOST_AUTO_TEST_SUITE(ON_GRID)

struct ON_GRID_PARAMETERS : public PARAMETERS {
  ON_GRID_PARAMETERS()
      : PARAMETERS(
            5,
            std::make_shared<DotVector>(DotVector{
                QuantumDot(Eigen::Vector3d::Zero(), Eigen::Vector3d(0, 0, 1)),
                QuantumDot(spacing * (num_boxes).cast<double>(),
                           Eigen::Vector3d(0, 0, 1))}))
  {
    history->fill(Eigen::Vector2cd::Zero());
    for(int i = -10; i < num_steps; ++i) {
      history->array[0][i][0](RHO_01) = src(i * dt);
      history->array[1][i][0](RHO_01) = 1;
    }
  }
};

BOOST_FIXTURE_TEST_CASE(GAUSSIAN_POINT_PROPAGATION, ON_GRID_PARAMETERS)
{
  AIM::AimInteraction aim(
      dots, history , interpolation_order, c, dt, grid, expansions,
      AIM::Expansions::Retardation(grid.max_transit_steps(c, dt) +
                                   interpolation_order),
      AIM::normalization::unit);

  const double delay =
      (dots->at(1).position() - dots->at(0).position()).norm() / c;

  const double err = test_propagation(
      aim, [&](const double t) { return src(t - delay); }, 1e-5);
  BOOST_TEST_MESSAGE("Maximum relative on-grid error: " << err);
}

BOOST_FIXTURE_TEST_CASE(GAUSSIAN_TIME_DERIVATIVE_PROPAGATION,
                        ON_GRID_PARAMETERS)
{
  AIM::AimInteraction aim(
      dots, history, interpolation_order, c, dt, grid, expansions,
      AIM::Expansions::TimeDerivative(grid.max_transit_steps(c, dt) +
                                      interpolation_order),
      AIM::normalization::unit);

  const double delay = 
      (dots->at(1).position() - dots->at(0).position()).norm() / c;

  double err = test_propagation(
      aim, [&](const double t) { return dt_src(t - delay); }, 1e-3);
  BOOST_TEST_MESSAGE("Maximum relative on-grid error: " << err);
}

BOOST_FIXTURE_TEST_CASE(GAUSSIAN_EFIE_PROPAGATION, ON_GRID_PARAMETERS)
{
  AIM::AimInteraction aim(
      dots, history, interpolation_order, c, dt, grid, expansions,
      AIM::Expansions::EFIE(grid.max_transit_steps(c, dt) + interpolation_order,
                            c),
      AIM::normalization::unit);

  using namespace std::placeholders;
  const auto dr = dots->at(1).position() - dots->at(0).position();
  double err = test_propagation(
      aim, std::bind(&ON_GRID_PARAMETERS::efie_src, this, dr, _1), 8e-1);
  BOOST_TEST_MESSAGE("Maximum relative on-grid error: " << err);
}

BOOST_AUTO_TEST_SUITE_END()  // ON_GRID

BOOST_AUTO_TEST_SUITE(OFF_GRID)

struct OFF_GRID_PARAMETERS : public PARAMETERS {
  OFF_GRID_PARAMETERS()
      : PARAMETERS(
            2,
            std::make_shared<DotVector>(DotVector{
                QuantumDot(spacing * Eigen::Vector3d(0.5, 0.5, 0.5).array(),
                           Eigen::Vector3d(0, 0, 1)),
                QuantumDot(spacing * ((num_boxes).cast<double>() + 0.5),
                           Eigen::Vector3d(0, 0, 1))}))
  {
    history->fill(Eigen::Vector2cd::Zero());
    for(int i = -10; i < num_steps; ++i) {
      history->array[0][i][0](RHO_01) = src(i * dt);
      history->array[1][i][0](RHO_01) = 1;
    }
  }
};

BOOST_FIXTURE_TEST_CASE(GAUSSIAN_POINT_PROPAGATION, OFF_GRID_PARAMETERS)
{
  AIM::AimInteraction aim(
      dots, history, interpolation_order, c, dt, grid, expansions,
      AIM::Expansions::Retardation(grid.max_transit_steps(c, dt) +
                                   interpolation_order),
      AIM::normalization::unit);

  const double delay = 
      (dots->at(1).position() - dots->at(0).position()).norm() / c;

  double err = test_propagation(
      aim, [&](const double t) { return src(t - delay); }, 1e-3);
  BOOST_TEST_MESSAGE("Maximum relative off-grid error: " << err);
}

BOOST_FIXTURE_TEST_CASE(GAUSSIAN_TIME_DERIVATIVE_PROPAGATION,
                        OFF_GRID_PARAMETERS)
{
  AIM::AimInteraction aim(
      dots, history, interpolation_order, c, dt, grid, expansions,
      AIM::Expansions::TimeDerivative(grid.max_transit_steps(c, dt) +
                                      interpolation_order),
      AIM::normalization::unit);

  const double delay = 
      (dots->at(1).position() - dots->at(0).position()).norm() / c;

  double err = test_propagation(
      aim, [&](const double t) { return dt_src(t - delay); }, 1e-1);
  BOOST_TEST_MESSAGE("Maximum relative off-grid error: " << err);
}

BOOST_FIXTURE_TEST_CASE(GAUSSIAN_EFIE_PROPAGATION, OFF_GRID_PARAMETERS)
{
  AIM::AimInteraction aim(
      dots, history, interpolation_order, c, dt, grid, expansions,
      AIM::Expansions::EFIE(grid.max_transit_steps(c, dt) + interpolation_order,
                            c),
      AIM::normalization::unit);

  using namespace std::placeholders;
  const auto dr = dots->at(1).position() - dots->at(0).position();
  double err = test_propagation(
      aim, std::bind(&OFF_GRID_PARAMETERS::efie_src, this, dr, _1), 20);
  BOOST_TEST_MESSAGE("Maximum relative off-grid error: " << err);
}

BOOST_AUTO_TEST_SUITE_END()  // OFF_GRID

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
  AIM::AimInteraction aim(dots, history, interp_order, c0, dt, grid,
                          expansions, nullptr, AIM::normalization::unit);
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
