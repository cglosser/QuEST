#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

#include "configuration.h"
#include "history.h"
#include "integrator.h"
#include "interactions/history_interaction.h"
#include "interactions/pulse_interaction.h"
#include "interactions/rotating_green_function.h"

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
    auto history(History::make_shared_history(config.num_particles, 22, config.num_timesteps));
    for(int t = -22; t <= 0 ; ++t) {
      for(int sol_idx = 0; sol_idx < config.num_particles; ++sol_idx) {
        (*history)[sol_idx][t][0] = Eigen::Vector2cd(1, 0); // Ground state
      }
    }

    auto pulse1 = make_shared<Pulse>(read_pulse_config(config.pulse_path));

    auto rotating_dyadic(make_shared<GreenFunction::RotatingDyadic>(
        config.mu0, config.c0, config.hbar, config.omega));

    auto hits = make_shared<HistoryInteraction>(qds, history, rotating_dyadic, config.interpolation_order);

    std::vector<std::shared_ptr<Interaction>> interactions{
        static_pointer_cast<Interaction>(make_shared<PulseInteraction>(qds, pulse1)),
        static_pointer_cast<Interaction>(hits)
    };

    PredictorCorrector::Integrator integrator(config.dt, 18, 22, 3.15, history,
                                              rhs_funs, std::move(interactions));

    cout << "Solving..." << endl;
    integrator.solve();

    cout << "Writing output..." << endl;

    ofstream output("output.dat");
    output << setprecision(12) << scientific;
    for(int t = 0; t < config.num_timesteps; ++t) {
      const double time = t * config.dt;
      output << time << " ";
      for(int sol_idx = 0; sol_idx < config.num_particles; ++sol_idx) {
        output << (*history)[sol_idx][t][0].transpose() << " ";
      }
      output << endl;
    }
  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
