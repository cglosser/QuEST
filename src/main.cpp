#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

#include "configuration.h"
#include "integrator.h"
#include "interaction_table.h"
#include "lagrange_set.h"
#include "math_utils.h"
#include "pulse.h"
#include "quantum_dot.h"

using namespace std;

int main(int argc, char *argv[])
{
  try {
    cout << setprecision(12) << scientific;
    auto vm = parse_configs(argc, argv);

    shared_ptr<DotVector> qds(new vector<QuantumDot>(2));
    (*qds)[0] = QuantumDot(Eigen::Vector3d(0, 0, 0), 2278.9013,
                           std::pair<double, double>(10, 20),
                           Eigen::Vector3d(5.2917721e-4, 0, 0));
    (*qds)[1] = QuantumDot(Eigen::Vector3d(0, 0, 0.010), 2278.9013,
                           std::pair<double, double>(10, 20),
                           Eigen::Vector3d(5.2917721e-4, 0, 0));

    const int num_steps = std::ceil(config.simulation_time / config.dt);
    auto superops(rhs_functions(*qds));

    shared_ptr<Pulse> p(new Pulse(15589.226 / 10, 5, 2278.9013, 2278.9013,
                                  Eigen::Vector3d(0, 0, 1),
                                  Eigen::Vector3d(1, 0, 0)));
    shared_ptr<Pulse> q(p);

    InteractionTable interaction_table(config.interpolation_order, qds, q);
    PredictorCorrector::Integrator solver(2, num_steps, config.dt, 18, 22, 3.15,
                                          superops, interaction_table);

    for(int t = -22; t <= num_steps; ++t) {
      solver.history[0][t][0] =
          Eigen::Vector2cd(0, gaussian((t * config.dt - 5)));
      solver.history[0][t][1] = Eigen::Vector2cd(0,0);

      if(t > 0) {
        solver.interaction_table.compute_history_interaction(solver.history, t);
      }

      cout << t * config.dt << " " << solver.history[0][t][0][1].real() << " "
           << solver.interaction_table.history_interaction[1] << endl;
    }

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
