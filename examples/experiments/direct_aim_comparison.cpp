#include <iomanip>
#include <iostream>

#include "interactions/AIM/aim_interaction.h"
#include "interactions/direct_interaction.h"
#include "math_utils.h"

struct PARAMETERS {
  int interpolation_order, expansion_order, num_steps, num_dots;
  double c, dt, total_time;

  Eigen::Array3i num_boxes;
  Eigen::Array3d spacing;

  std::shared_ptr<DotVector> dots;
  std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history;

  AIM::Grid grid;
  AIM::Expansions::ExpansionTable expansions;

  PARAMETERS()
      : interpolation_order(4),
        expansion_order(1),
        num_steps(1024),
        num_dots(2),

        c(1),
        dt(1),
        total_time(num_steps * dt),

        num_boxes(Eigen::Vector3i(8, 8, 8)),
        spacing(Eigen::Array3d(1, 1, 1) * c * dt),

        dots(std::make_shared<DotVector>(DotVector{
            QuantumDot(Eigen::Vector3d::Zero(), Eigen::Vector3d(0, 0, 1)),
            QuantumDot(spacing * num_boxes.array().cast<double>(),
                       Eigen::Vector3d(0, 0, 1))})),
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
      Propagation::RotatingFramePropagator(1, params.c, 1, 20 * M_PI),
      params.interpolation_order, params.c, params.dt);

  for(int i = 0; i < params.num_steps; ++i) {
    std::cout << i << " " << direct.evaluate(i).transpose() << std::endl;
  }

  return 0;
}
