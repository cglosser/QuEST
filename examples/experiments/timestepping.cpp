#include <iomanip>
#include <iostream>
#include <limits>

#include "interactions/AIM/aim_interaction.h"
#include "interactions/direct_interaction.h"
#include "math_utils.h"

using dbl = std::numeric_limits<double>;

struct PARAMETERS {
  int interpolation_order, expansion_order, num_dots;
  double c, omega;

  double total_time, dt;
  int num_steps;

  Eigen::Array3i num_boxes;
  Eigen::Array3d spacing;

  std::shared_ptr<DotVector> dots;
  std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history;

  AIM::Grid grid;
  AIM::Expansions::ExpansionTable expansions;

  PARAMETERS()
      : interpolation_order(4),
        expansion_order(4),
        num_dots(2),

        c(1),
        omega(M_PI / 10),

        total_time(1024),
        dt(1),
        num_steps(total_time / dt + 1),

        num_boxes(Eigen::Vector3i(8, 8, 8)),
        spacing(Eigen::Array3d(1, 1, 1) * c * M_PI / (10 * omega)),

        dots(std::make_shared<DotVector>(DotVector{
            QuantumDot(Eigen::Vector3d(0, 0, 0), Eigen::Vector3d(0, 0, 1)),
            QuantumDot(Eigen::Vector3d(10.5, 10.5, 10.5), Eigen::Vector3d(0, 0, 1))})),
        history(std::make_shared<Integrator::History<Eigen::Vector2cd>>(
            num_dots, 10, num_steps)),

        grid(spacing, dots, expansion_order),
        expansions(AIM::Expansions::LeastSquaresExpansionSolver::get_expansions(
            expansion_order, grid, *dots))
  {
    history->fill(Eigen::Vector2cd::Zero());
    for(int i = -10; i < num_steps; ++i) {
      history->array[0][i][0](RHO_01) = gaussian_source(i * dt);
    }
  };

  double gaussian_source(const double t) const
  {
    double arg = (t - total_time / 2.0) / (total_time / 12.0);
    return gaussian(arg);
  }
};

int main()
{
  PARAMETERS params;

  DirectInteraction direct(
      params.dots, params.history,
      Propagation::RotatingFramePropagator(1, params.c, params.omega),
      params.interpolation_order, params.c, params.dt);

  AIM::Grid grid(params.spacing, params.dots, params.expansion_order);
  auto expansion_table =
      AIM::Expansions::LeastSquaresExpansionSolver::get_expansions(
          params.expansion_order, grid, *params.dots);

  AIM::AimInteraction aim(
      params.dots, params.history, params.interpolation_order, params.c,
      params.dt, grid, expansion_table,
      AIM::Expansions::RotatingEFIE(
          grid.max_transit_steps(params.c, params.dt) +
              params.interpolation_order,
          params.c, params.dt, params.omega),
      AIM::normalization::Helmholtz(params.omega / params.c, 1));

  std::cout.precision(dbl::max_digits10);
  for(int i = 0; i < params.num_steps; ++i) {
    std::cout << i * params.dt << " " << direct.evaluate(i)(params.num_dots - 1) << " "
              << aim.evaluate(i)(params.num_dots - 1) << std::endl;
  }

  return 0;
}
