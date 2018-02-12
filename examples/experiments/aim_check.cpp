#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "interactions/AIM/aim_interaction.h"
#include "interactions/direct_interaction.h"
#include "math_utils.h"

using dbl = std::numeric_limits<double>;

class Gaussian {
 public:
  Gaussian(double mu, double sigma) : mu(mu), sigma(sigma){};
  double operator()(double t) const { return gaussian((t - mu) / sigma); }
 private:
  double mu, sigma;
};

int main()
{
  const int interpolation_order = 4, expansion_order = 3, num_steps = 2048;
  const double c = 1, dt = 0.5, total_time = num_steps * dt, omega = 0;
  const Gaussian source(total_time / 2.0, total_time / 12.0);

  const Eigen::Array3d spacing(1.6, 1.6, 1.6);

  auto dots = std::make_shared<DotVector>();
  dots->push_back(QuantumDot({0.0, 0.0, 0.0}, omega, {0, 0}, {0, 0, 1}));
  dots->push_back(QuantumDot({0.0, 0.0, 10.0}, omega, {0, 0}, {0, 0, 1}));

  const int num_dots = dots->size();
  auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(
      num_dots, 10, num_steps);
  for(int d = 0; d < num_dots; ++d) {
    for(int t = -10; t < num_steps; ++t) {
      history->array_[d][t][0](RHO_01) = source(t * dt);
    }
  }

  Propagation::Laplace<cmplx> laplace;

  DirectInteraction direct(dots, history, laplace, interpolation_order, c, dt);

  std::ofstream direct_fd("direct.dat");
  direct_fd.precision(17);
  for(int t = 0; t < num_steps; ++t) {
    direct_fd << direct.evaluate(t).transpose() << std::endl;
  }

  AIM::Grid grid(spacing, expansion_order, *dots);
  auto expansions =
      AIM::Expansions::LeastSquaresExpansionSolver::get_expansions(
          expansion_order, grid, *dots);

  AIM::DirectInteraction di(dots, history, laplace, interpolation_order, 1, c,
                            dt, grid);
  AIM::Nearfield nf(dots, history, interpolation_order, 1, c, dt, grid,
                    expansions,
                    AIM::Expansions::Retardation(grid.max_transit_steps(c, dt) +
                                                 interpolation_order),
                    AIM::normalization::Laplace());
  AIM::Farfield ff(dots, history, interpolation_order, c, dt, grid, expansions,
                   AIM::Expansions::Retardation(grid.max_transit_steps(c, dt) +
                                                interpolation_order),
                   AIM::normalization::Laplace());

  std::ofstream comparison("comparison.dat");
  comparison.precision(17);
  for(int t = 0; t < num_steps; ++t) {
    std::cout << t << std::endl;
    comparison << t * dt << " "
               << (direct.evaluate(t) -
                   (di.evaluate(t) + (ff.evaluate(t) - nf.evaluate(t))))
                      .matrix()
                      .norm()
               << std::endl;
  }

  return 0;
}
