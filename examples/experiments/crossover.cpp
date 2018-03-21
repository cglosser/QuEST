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

const std::vector<Eigen::Vector3d> pos8 = {
    {0.25, 0.25, 0.25}, {0.25, 0.25, 0.75}, {0.25, 0.75, 0.25},
    {0.25, 0.75, 0.75}, {0.75, 0.25, 0.25}, {0.75, 0.25, 0.75},
    {0.75, 0.75, 0.25}, {0.75, 0.75, 0.75}};

const std::vector<Eigen::Vector3d> pos27 = {
    {0.166667, 0.166667, 0.166667}, {0.166667, 0.166667, 0.5},
    {0.166667, 0.166667, 0.833333}, {0.166667, 0.5, 0.166667},
    {0.166667, 0.5, 0.5},           {0.166667, 0.5, 0.833333},
    {0.166667, 0.833333, 0.166667}, {0.166667, 0.833333, 0.5},
    {0.166667, 0.833333, 0.833333}, {0.5, 0.166667, 0.166667},
    {0.5, 0.166667, 0.5},           {0.5, 0.166667, 0.833333},
    {0.5, 0.5, 0.166667},           {0.5, 0.5, 0.5},
    {0.5, 0.5, 0.833333},           {0.5, 0.833333, 0.166667},
    {0.5, 0.833333, 0.5},           {0.5, 0.833333, 0.833333},
    {0.833333, 0.166667, 0.166667}, {0.833333, 0.166667, 0.5},
    {0.833333, 0.166667, 0.833333}, {0.833333, 0.5, 0.166667},
    {0.833333, 0.5, 0.5},           {0.833333, 0.5, 0.833333},
    {0.833333, 0.833333, 0.166667}, {0.833333, 0.833333, 0.5},
    {0.833333, 0.833333, 0.833333}};

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
const Eigen::Vector3d spacing(c *dt, c *dt, c *dt);
const int interpolation_order = 5, expansion_order = 1;

DotVector make_system(const int num_boxes)
{
  DotVector dots;
  for(int x = 0; x < num_boxes + 1; ++x) {
    for(int y = 0; y < num_boxes + 1; ++y) {
      for(int z = 0; z < num_boxes + 1; ++z) {
        for(const auto &r : pos8) {
          dots.push_back(QuantumDot(r + Eigen::Vector3d(x, y, z), {0, 0, 1}));
        }
      }
    }
  }

  return dots;
}

int main()
{
  using hist_t = Integrator::History<Eigen::Vector2cd>;
  using interval_t = std::chrono::duration<double>;

  std::ofstream fd("timing.dat");
  fd << std::setprecision(17);

  for(int n_boxes = 0; n_boxes < 6; ++n_boxes) {
    const auto dots = std::make_shared<DotVector>(make_system(n_boxes));
    auto history = std::make_shared<hist_t>(dots->size(), 10, num_steps);
    for(int dot = 0; dot < 8; ++dot) {
      for(int t = -10; t < num_steps; ++t) {
        history->array_[dot][t][0](RHO_01) =
            1e-8 * gaussian((t * dt - num_steps / 2.0) / (num_steps / 12.0));
      }
    }

    // == SLOW STUFF =================================================

    Propagation::Laplace<cmplx> kernel;
    DirectInteraction direct(dots, history, kernel, interpolation_order, c, dt);

    auto slow_start = std::chrono::high_resolution_clock::now();
    for(int t = 0; t < num_steps; ++t) {
      direct.evaluate(t);
    }
    auto slow_stop = std::chrono::high_resolution_clock::now();
    interval_t slow_span =
        std::chrono::duration_cast<interval_t>(slow_stop - slow_start);

    // == FAST STUFF =================================================

    // FFT
    std::ofstream far("farfield.dat");

    using LSE = AIM::Expansions::LeastSquaresExpansionSolver;

    auto grid = std::make_shared<AIM::Grid>(spacing, expansion_order, *dots);
    auto expansion_table = std::make_shared<AIM::Expansions::ExpansionTable>(
        LSE::get_expansions(expansion_order, *grid, *dots));

    AIM::Farfield ff(dots, history, interpolation_order, c, dt, grid,
                     expansion_table,
                     AIM::Expansions::Retardation(
                         grid->max_transit_steps(c, dt) + interpolation_order),
                     AIM::Normalization::Laplace());

    auto fast_start = std::chrono::high_resolution_clock::now();
    for(int t = 0; t < num_steps; ++t) {
      ff.evaluate(t);
    }
    auto fast_stop = std::chrono::high_resolution_clock::now();
    interval_t fast_span =
        std::chrono::duration_cast<interval_t>(fast_stop - fast_start);

    // Nearfield

    AIM::Nearfield nf(dots, history, interpolation_order, c, dt, grid,
                      expansion_table, nullptr, AIM::Normalization::Laplace(),
                      std::make_shared<std::vector<AIM::Grid::ipair_t>>(
                          grid->nearfield_point_pairs(1, *dots)),
                      0);

    auto near_start = std::chrono::high_resolution_clock::now();
    for(int t = 0; t < num_steps; ++t) {
      nf.evaluate(t);
    }
    auto near_stop = std::chrono::high_resolution_clock::now();
    interval_t near_span =
        std::chrono::duration_cast<interval_t>(near_stop - near_start);

    // == OUTPUT =====================================================

    std::cout << n_boxes << " " << dots->size() << " " << slow_span.count()
              << " " << near_span.count() << " " << fast_span.count()
              << std::endl;
    fd << n_boxes << " " << dots->size() << " " << slow_span.count() << " "
       << near_span.count() << " " << fast_span.count() << std::endl;
  }
  return 0;
}
