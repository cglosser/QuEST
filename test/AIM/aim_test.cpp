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

BOOST_AUTO_TEST_SUITE_END()  // AIM
