#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

#include "configuration.h"
#include "integrator.h"
#include "interactions/pulse_interaction.h"
#include "pulse.h"

using namespace std;

int main(int argc, char *argv[])
{
  try {
    cout << setprecision(12) << scientific;
    auto vm = parse_configs(argc, argv);

    auto qds(make_shared<DotVector>(config.num_particles));
    (*qds)[0] = QuantumDot(Eigen::Vector3d(0, 0, 0), 2278.9013,
                           std::pair<double, double>(10, 10),
                           Eigen::Vector3d(5.2917721e-4, 0, 0));
    auto rhs_funs(rhs_functions(*qds, config.omega));

    // Set up History
    typedef History::HistoryArray::extent_range range;
    auto history(make_shared<History::HistoryArray>(boost::extents[config.num_particles][range(-22, 1000)][2]));
    for(int t = -22; t <= 0 ; ++t) {
      (*history)[0][t][0] = Eigen::Vector2cd(1, 0); // Ground state
    }

    // Set up Pulse and Interaction
    auto pulse(
        make_shared<Pulse>(1558.92260227, 5, config.omega, config.omega,
                           Eigen::Vector3d(0, 0, config.omega / config.c0),
                           Eigen::Vector3d(1, 0, 0)));

    auto field_interarction(make_shared<PulseInteraction>(qds, pulse));

    // Set up Integrator
    PredictorCorrector::Integrator integrator(config.dt, 18, 22, 3.15, history,
                                              rhs_funs, field_interarction);

    // Solve the system
    for(int i = 0; i < 1000; ++i) {
      integrator.solve(i);
    }

    for(int t = 0; t < 1000; ++t) {
      const double time = t * config.dt;
      cout << time << " " << (*history)[0][t][0].transpose() << endl;
    }
  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
