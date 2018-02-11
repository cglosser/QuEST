#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

#include "configuration.h"
#include "integrator/RHS/bloch_rhs.h"
#include "integrator/history.h"
#include "integrator/integrator.h"
#include "interactions/direct_interaction.h"
#include "interactions/green_function.h"
#include "interactions/pulse_interaction.h"

using namespace std;

int main(int argc, char *argv[])
{
  try {
    cout << setprecision(12) << scientific;
    auto vm = parse_configs(argc, argv);

    cout << "Initializing..." << endl;

    auto qds = make_shared<DotVector>(import_dots(config.qd_path));
    qds->resize(config.num_particles);
    auto rhs_funs = rhs_functions(*qds, config.omega);

    // Set up History
    auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(
        config.num_particles, 22, config.num_timesteps);
    history->fill(Eigen::Vector2cd(0, 0));
    history->initialize_past(Eigen::Vector2cd(1, 0));

    // Set up Interactions
    auto pulse1 = make_shared<Pulse>(read_pulse_config(config.pulse_path));
    Propagation::RotatingEFIE rotating_dyadic(
        config.c0, config.mu0 / (4 * M_PI * config.hbar), config.omega);

    std::vector<std::shared_ptr<InteractionBase>> interactions{
        make_shared<PulseInteraction>(qds, pulse1, config.hbar, config.dt),
        make_shared<DirectInteraction>(qds, history, rotating_dyadic,
                                       config.interpolation_order, config.dt,
                                       config.c0)};

    // Set up RHS functions
    auto rhs_funcs = rhs_functions(*qds, config.omega);

    // Set up Bloch RHS
    std::unique_ptr<Integrator::RHS<Eigen::Vector2cd>> bloch_rhs =
        std::make_unique<Integrator::BlochRHS>(
            config.dt, history, std::move(interactions), std::move(rhs_funcs));

    Integrator::PredictorCorrector<Eigen::Vector2cd> solver(
        config.dt, 18, 22, 3.15, history, std::move(bloch_rhs));

    cout << "Solving..." << endl;
    solver.solve();

    cout << "Writing output..." << endl;
    ofstream outfile("output.dat");
    outfile << scientific << setprecision(15);
    for(int t = 0; t < config.num_timesteps; ++t) {
      for(int n = 0; n < config.num_particles; ++n) {
        outfile << history->array_[n][t][0].transpose() << " ";
      }
      outfile << "\n";
    }

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
