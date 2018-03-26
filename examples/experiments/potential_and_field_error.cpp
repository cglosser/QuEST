#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "interactions/AIM/aim_interaction.h"
#include "interactions/direct_interaction.h"
#include "math_utils.h"

using dbl = std::numeric_limits<double>;

class Gaussian {
 public:
  Gaussian(double mu, double sigma) : mu_(mu), sigma_(sigma){};
  double operator()(double t) const { return gaussian((t - mu_) / sigma_); }
 private:
  double mu_, sigma_;
};

const std::vector<Eigen::Vector3d> pos64 = {
    {0.125, 0.125, 0.125}, {0.125, 0.125, 0.375}, {0.125, 0.125, 0.625},
    {0.125, 0.125, 0.875}, {0.125, 0.375, 0.125}, {0.125, 0.375, 0.375},
    {0.125, 0.375, 0.625}, {0.125, 0.375, 0.875}, {0.125, 0.625, 0.125},
    {0.125, 0.625, 0.375}, {0.125, 0.625, 0.625}, {0.125, 0.625, 0.875},
    {0.125, 0.875, 0.125}, {0.125, 0.875, 0.375}, {0.125, 0.875, 0.625},
    {0.125, 0.875, 0.875}, {0.375, 0.125, 0.125}, {0.375, 0.125, 0.375},
    {0.375, 0.125, 0.625}, {0.375, 0.125, 0.875}, {0.375, 0.375, 0.125},
    {0.375, 0.375, 0.375}, {0.375, 0.375, 0.625}, {0.375, 0.375, 0.875},
    {0.375, 0.625, 0.125}, {0.375, 0.625, 0.375}, {0.375, 0.625, 0.625},
    {0.375, 0.625, 0.875}, {0.375, 0.875, 0.125}, {0.375, 0.875, 0.375},
    {0.375, 0.875, 0.625}, {0.375, 0.875, 0.875}, {0.625, 0.125, 0.125},
    {0.625, 0.125, 0.375}, {0.625, 0.125, 0.625}, {0.625, 0.125, 0.875},
    {0.625, 0.375, 0.125}, {0.625, 0.375, 0.375}, {0.625, 0.375, 0.625},
    {0.625, 0.375, 0.875}, {0.625, 0.625, 0.125}, {0.625, 0.625, 0.375},
    {0.625, 0.625, 0.625}, {0.625, 0.625, 0.875}, {0.625, 0.875, 0.125},
    {0.625, 0.875, 0.375}, {0.625, 0.875, 0.625}, {0.625, 0.875, 0.875},
    {0.875, 0.125, 0.125}, {0.875, 0.125, 0.375}, {0.875, 0.125, 0.625},
    {0.875, 0.125, 0.875}, {0.875, 0.375, 0.125}, {0.875, 0.375, 0.375},
    {0.875, 0.375, 0.625}, {0.875, 0.375, 0.875}, {0.875, 0.625, 0.125},
    {0.875, 0.625, 0.375}, {0.875, 0.625, 0.625}, {0.875, 0.625, 0.875},
    {0.875, 0.875, 0.125}, {0.875, 0.875, 0.375}, {0.875, 0.875, 0.625},
    {0.875, 0.875, 0.875}};

const int num_steps = 1024;
const double c = 1, omega = 0;
const double dt = 1, total_time = dt * num_steps;
const Eigen::Vector3d spacing = Eigen::Vector3d(c * dt, c *dt, c *dt) / 2;
const int interpolation_order = 5, expansion_order = 5;

DotVector make_system(const Eigen::Vector3d &dr)
{
  DotVector dots;

  for(const auto &p : pos64) {
    dots.emplace_back(p, Eigen::Vector3d(0, 0, 1));
  }

  for(const auto &p : pos64) {
    dots.emplace_back(p + dr, Eigen::Vector3d(0, 0, 1));
  }

  return dots;
}

int main()
{
  using hist_t = Integrator::History<Eigen::Vector2cd>;

  const auto dots = std::make_shared<DotVector>(make_system({10, 10, 10}));

  const Gaussian source(total_time / 2.0, total_time / 12.0);

  const int num_dots = dots->size();
  auto history = std::make_shared<hist_t>(num_dots, 10, num_steps);
  for(int t = -10; t < num_steps; ++t) {
    history->array_[0][t][0](RHO_01) = source(t * dt);
  }

  // == Direct =======================================================

  Propagation::DelSq_Laplace<cmplx> kernel(c);
  DirectInteraction direct(dots, history, kernel, interpolation_order, c, dt);

  // == AIM ==========================================================

  using LSE = AIM::Expansions::LeastSquaresExpansionSolver;

  auto grid = std::make_shared<AIM::Grid>(spacing, expansion_order, *dots);
  auto expansion_table = std::make_shared<AIM::Expansions::ExpansionTable>(
      LSE::get_expansions(expansion_order, *grid, *dots));

  std::cout << "Shape: " << grid->shape().transpose() << std::endl;

  AIM::Expansions::Del_Del exp_fun(grid->max_transit_steps(c, dt) +
                                   interpolation_order);

  AIM::Normalization::Laplace norm_fun;

  AIM::Farfield ff(dots, history, interpolation_order, c, dt, grid,
                   expansion_table, exp_fun, norm_fun);

  // == Evolution ====================================================

  std::array<std::ofstream, 2> fd{std::ofstream("aim.dat"),
                                  std::ofstream("direct.dat")};
  for(auto &f : fd) f.precision(17);

  for(int t = 0; t < num_steps; ++t) {
    fd[0] << ff.evaluate(t).transpose().real() << std::endl;
    fd[1] << direct.evaluate(t).transpose().real() << std::endl;
  }

  return 0;
}
