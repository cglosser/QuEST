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
  double operator()(double t) const { return Math::gaussian((t - mu) / sigma); }
 private:
  double mu, sigma;
};

int main()
{
  const int interpolation_order = 4, expansion_order = 1, num_steps = 1024;
  const double c = 299.792458, dt = 4.136e-2, total_time = dt * num_steps,
               omega = 2278.9013;
  const Gaussian source(total_time / 2.0, total_time / 12.0);

  const double wavelength = 2 * M_PI * c / omega;

  const double h = wavelength / 5;
  const Eigen::Array3d spacing = Eigen::Array3d(h, h, h);

  const Eigen::Array3d pos(wavelength / 10, wavelength / 10, wavelength / 10);

  auto dots = std::make_shared<DotVector>();
  dots->push_back(QuantumDot(pos, 2278.9013, {10, 20}, {0.00052917721, 0, 0}));
  dots->push_back(QuantumDot(pos + Eigen::Array3d(0, 0, 2 * wavelength),
                             2278.9013, {10, 20}, {0.00052917721, 0, 0}));

  const int num_dots = dots->size();
  auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(
      num_dots, 10, num_steps);
  for(int t = -10; t < num_steps; ++t) {
    history->array_[0][t][0](RHO_01) = source(t * dt);
  }

  Propagation::Laplace<cmplx> prop(2.4341349e-5);
  auto grid = AIM::Grid(spacing, expansion_order, *dots);
  int transit = grid.max_transit_steps(c, dt) + interpolation_order;

  AIM::Interaction aim(dots, history, prop, spacing, interpolation_order,
                       expansion_order, 1, c, dt,
                       AIM::Normalization::Laplace(2.4341349e-5),
                       Projector::Potential<cmplx>(transit));
  DirectInteraction direct(dots, history, prop, interpolation_order, c, dt);

  std::ofstream aim_fd("aim.dat"), dir_fd("direct.dat");
  aim_fd.precision(17);
  dir_fd.precision(17);

  std::cout << "Shape: " << grid.shape().transpose() << std::endl;
  for(int t = 0; t < num_steps; ++t) {
    if(t % 100 == 0) std::cout << t << std::endl;
    auto dir = direct.evaluate(t);
    dir_fd << dir.transpose() << "\n";

    auto aux = aim.evaluate(t);
    aim_fd << aux.transpose() << "\n";
  }

  return 0;
}
