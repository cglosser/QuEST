#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "interactions/AIM/aim_interaction.h"
#include "interactions/direct_interaction.h"
#include "lagrange_set.h"
#include "math_utils.h"

using dbl = std::numeric_limits<double>;

class Gaussian {
 public:
  Gaussian(double mu, double sigma) : mu_(mu), sigma_(sigma){};
  double operator()(double t) const { return gaussian((t - mu_) / sigma_); }
 private:
  double mu_, sigma_;
};

class Logistic {
 public:
  Logistic(double mu, double sigma) : mu_(mu), sigma_(sigma){};
  double operator()(double t) const
  {
    return 1 / (1 + std::exp(-(t - mu_) / sigma_));
  }

 private:
  double mu_, sigma_;
};

int main()
{
  const int num_steps = 1024;

  const double c = 299.792458, omega = 2278.9013, k2 = 2.4341348e-05;
  const double dt = 0.5e-2, total_time = dt * num_steps;

  // const double c = 1, omega = 0, k2 = 1;
  // const double dt = 1, total_time = dt * num_steps;

  const int interpolation_order = 5, expansion_order = 6;

  const double s = c * 2 * M_PI / (10 * omega);
  const Eigen::Array3d spacing(s, s, s);

  auto dots = std::make_shared<DotVector>();
  dots->push_back(QuantumDot(Eigen::Vector3d(0, 0, 0), {0, 1, 0}));
  dots->push_back(QuantumDot(Eigen::Vector3d(0, 0, 0.867825), {0, 1, 0}));

  const Gaussian source(total_time / 2.0, total_time / 12.0);

  const int num_dots = dots->size();
  auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(
      num_dots, 10, num_steps);
  for(int t = -10; t < num_steps; ++t) {
    history->array_[0][t][0](RHO_01) = source(t * dt);
  }
  // == Direct =======================================================

  Propagation::EFIE<cmplx> kernel(c, k2);
  DirectInteraction direct(dots, history, kernel, interpolation_order, c, dt);

  // == AIM ==========================================================

  using LSE = AIM::Expansions::LeastSquaresExpansionSolver;

  auto grid = std::make_shared<AIM::Grid>(spacing, expansion_order, *dots);
  auto expansion_table = std::make_shared<AIM::Expansions::ExpansionTable>(
      LSE::get_expansions(expansion_order, *grid, *dots));

  std::cout << "Shape: " << grid->shape().transpose() << std::endl;

  AIM::Expansions::RotatingEFIE exp_fun(
      grid->max_transit_steps(c, dt) + interpolation_order, c, dt, omega);

  AIM::Normalization::Laplace norm_fun(k2);

  AIM::Farfield ff(dots, history, interpolation_order, c, dt, grid,
                   expansion_table, exp_fun, norm_fun);

  AIM::Nearfield nf(dots, history, interpolation_order, c, dt, grid,
                    expansion_table, nullptr, norm_fun,
                    std::make_shared<std::vector<AIM::Grid::ipair_t>>(
                        grid->nearfield_point_pairs(100, *dots)),
                    omega);

  AIM::Interaction combined(dots, history, kernel, spacing, interpolation_order,
                            expansion_order, 1, c, dt, exp_fun, norm_fun,
                            omega);

  // =================================================================

  std::array<std::ofstream, 4> fd{
      std::ofstream("nearfield.dat"), std::ofstream("farfield.dat"),
      std::ofstream("direct.dat"), std::ofstream("combined.dat")};
  for(auto &f : fd) f.precision(17);

  for(int t = 0; t < num_steps; ++t) {
    if(t % 100 == 0) std::cout << t << std::endl;
    fd[0] << nf.evaluate(t).transpose() << std::endl;
    fd[1] << ff.evaluate(t).conjugate().transpose() << std::endl;
    fd[2] << direct.evaluate(t).transpose() << std::endl;
    fd[3] << combined.evaluate(t).transpose() << std::endl;
  }

  return 0;
}
