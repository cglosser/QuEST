// This file looks an _awful_ lot like the unit tests for testing the
// propagation + expansion functionality! That's probably because a lot of this 
// came directly from there. Be warned, however: that file is the "canonical"
// expression of this functionality and has more rigorous checks in place to validate
// what's going on. Don't use this file unless those tests pass *and* you know
// what you're doing here.
//                                                              --Connor

#include <iomanip>
#include <iostream>

#include "interactions/AIM/aim_interaction.h"
#include "math_utils.h"

struct PARAMETERS {
  static int interpolation_order, num_steps, num_dots;
  static double c, dt, total_time;

  static Eigen::Array3i num_boxes;
  static Eigen::Array3d spacing;

  int expansion_order;  // Different orders for different test geometries

  std::shared_ptr<DotVector> dots;
  std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history;

  AIM::Grid grid;
  AIM::Expansions::ExpansionTable expansions;

  PARAMETERS(const int expansion_order, std::shared_ptr<DotVector> dots)
      : expansion_order(expansion_order),
        dots(dots),
        history(std::make_shared<Integrator::History<Eigen::Vector2cd>>(
            num_dots, 10, num_steps)),
        grid(spacing, dots, expansion_order),
        expansions(AIM::Expansions::LeastSquaresExpansionSolver::get_expansions(
            expansion_order, grid, *dots)){};

  virtual ~PARAMETERS() = 0;

  double gaussian_source(const double t) const
  {
    double arg = (t - total_time / 2.0) / (total_time / 12.0);
    return gaussian(arg);
  }

  double gaussian_time_derivative(const double t) const
  {
    const double mu = total_time / 2.0, sigma = total_time / 12.0;
    const double arg = (t - mu) / sigma;
    return -arg * gaussian(arg) / sigma;
  }

  double gaussian_efie(const Eigen::Vector3d &dr, const double t)
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

  void output_propagation_data(AIM::AimInteraction &aim,
                               const std::function<double(const double)> &sol,
                               const std::string &fname)
  {
    std::ofstream output(fname, std::ofstream::out);
    output << std::setprecision(17);
    output << std::scientific;
    output << "Time,Input,Output (analytic),Output (calculated)" << std::endl;

    for(int i = 0; i < num_steps; ++i) {
      auto x = aim.evaluate(i);
      output << i * dt << ", " << gaussian_source(i * dt) << ", " << sol(i * dt)
             << ", " << x(1).real() << std::endl;
    }
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
      history->array[0][i][0](RHO_01) = gaussian_source(i * dt);
      history->array[1][i][0](RHO_01) = 1;
    }
  }

  void retardation_test(const std::string &fname)
  {
    AIM::AimInteraction aim(
        dots, history, nullptr, interpolation_order, c, dt, grid, expansions,
        AIM::Expansions::Retardation(grid.max_transit_steps(c, dt) +
                                     interpolation_order),
        AIM::normalization::unit);

    const double delay =
        (dots->at(1).position() - dots->at(0).position()).norm() / c;

    output_propagation_data(
        aim, [&](const double t) { return gaussian_source(t - delay); }, fname);
  }

  void time_derivative_test(const std::string &fname)
  {
    AIM::AimInteraction aim(
        dots, history, nullptr, interpolation_order, c, dt, grid, expansions,
        AIM::Expansions::TimeDerivative(grid.max_transit_steps(c, dt) +
                                        interpolation_order),
        AIM::normalization::unit);

    const double delay =
        (dots->at(1).position() - dots->at(0).position()).norm() / c;

    output_propagation_data(
        aim,
        [&](const double t) { return gaussian_time_derivative(t - delay); },
        fname);
  }

  void efie_test(const std::string &fname)
  {
    AIM::AimInteraction aim(
        dots, history, nullptr, interpolation_order, c, dt, grid, expansions,
        AIM::Expansions::EFIE(
            grid.max_transit_steps(c, dt) + interpolation_order, c),
        AIM::normalization::unit);

    using namespace std::placeholders;
    const auto dr = dots->at(1).position() - dots->at(0).position();
    output_propagation_data(
        aim, std::bind(&ON_GRID_PARAMETERS::gaussian_efie, this, dr, _1),
        fname);
  }
};

int main()
{
  ON_GRID_PARAMETERS on_grid;
  on_grid.retardation_test("gaussian_retardation.dat");
  on_grid.time_derivative_test("gaussian_time_derivative.dat");
  on_grid.efie_test("gaussian_efie.dat");

  return 0;
}
