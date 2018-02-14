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
  const int interpolation_order = 4, expansion_order = 5, num_steps = 1024;
  const double c = 299.792458, dt = 4.136e-2, total_time = dt * num_steps,
               omega = 2278.9013;
  const Gaussian source(total_time / 2.0, total_time / 12.0);

  const double h = 2 * M_PI * c / (20 * omega);
  const Eigen::Array3d spacing(h, h, h);

  auto dots = std::make_shared<DotVector>();
  dots->push_back(
      QuantumDot({0, 0, 0}, 2278.9013, {10, 20}, {0.00052917721, 0, 0}));
  dots->push_back(QuantumDot({h / 2, h / 2, h / 2}, 2278.9013, {10, 20},
                             {0.00052917721, 0, 0}));

  // std::cout << omega << std::endl;
  // std::cout << h << std::endl;
  // std::cout << (dots->at(1).position() - dots->at(0).position()).norm()
  //<< std::endl;
  // std::cout << dt << std::endl;
  // std::cout << total_time << std::endl;

  const int num_dots = dots->size();
  auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(
      num_dots, 10, num_steps);
  for(int t = -10; t < num_steps; ++t) {
    history->array_[0][t][0](RHO_01) = source(t * dt);
  }

  Propagation::RotatingEFIE prop(c, 2.4341349e-5, omega);
  auto grid = AIM::Grid(spacing, expansion_order, *dots);
  int nt = grid.max_transit_steps(c, dt) + interpolation_order;

  AIM::Interaction aim(dots, history, prop, spacing, interpolation_order,
                       expansion_order, 1, c, dt,
                       AIM::Expansions::RotatingEFIE(nt, c, dt, omega),
                       AIM::normalization::Laplace(2.4341349e-5));
  DirectInteraction direct(dots, history, prop, interpolation_order, c, dt);

  std::ofstream comparison("comparison.dat");
  comparison.precision(17);
  for(int t = 0; t < num_steps; ++t) {
    if(t % 100 == 0) std::cout << t << std::endl;
    comparison << direct.evaluate(t).transpose() << " "
               << aim.evaluate(t).transpose() << "\n";
  }

  return 0;
}
