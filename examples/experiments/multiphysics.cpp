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
  const int interpolation_order = 4, expansion_order = 4, num_steps = 4000;
  const double c = 299.792458, dt = 0.5e-2, total_time = 20, omega = 2278.9013;
  const Gaussian source(total_time / 2.0, total_time / 12.0);

  const Eigen::Array3d spacing(0.10, 0.10, 0.10);

  auto dots = std::make_shared<DotVector>();
  dots->push_back(QuantumDot(
      {-0.006772025629369738, 0.006193253263435561, -0.010766531732233575},
      2278.9013, {10, 20}, {0.00052917721, 0, 0}));
  dots->push_back(QuantumDot(
      {0.011823480098434758, 0.010089453212510147, -0.008611542585243048},
      2278.9013, {10, 20}, {0.00052917721, 0, 0}));

  const int num_dots = dots->size();
  auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(
      num_dots, 10, num_steps);
  for(int d = 0; d < num_dots; ++d) {
    for(int t = -10; t < num_steps; ++t) {
      history->array_[d][t][0](RHO_01) = source(t * dt);
    }
  }

  Propagation::RotatingEFIE efie(c, 2.4341349e-5, omega);
  int nt = AIM::Grid(spacing, expansion_order, *dots).max_transit_steps(c, dt) +
           interpolation_order;

  AIM::Interaction aim(dots, history, efie, spacing, interpolation_order,
                       expansion_order, 1, c, dt,
                       AIM::Expansions::RotatingEFIE(nt, c, dt, omega),
                       AIM::normalization::Helmholtz(omega / c, 2.4341349e-5));
  DirectInteraction direct(dots, history, efie, interpolation_order, c, dt);

  std::ofstream comparison("comparison.dat");
  comparison.precision(17);
  for(int t = 0; t < num_steps; ++t) {
    std::cout << t << std::endl;
    // comparison << t * dt << " "
    //<< (direct.evaluate(t) - aim.evaluate(t)).matrix().norm()
    //<< std::endl;
    comparison << direct.evaluate(t).transpose() << "     " << aim.evaluate(t).transpose() << std::endl;
  }

  return 0;
}
