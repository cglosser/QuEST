#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

#include "configuration.h"
#include "interactions/history_interaction.h"
#include "interactions/rotating_green_function.h"
#include "math_utils.h"

using namespace std;

int main(int argc, char *argv[])
{
  try {
    cout << setprecision(12) << scientific;
    auto vm = parse_configs(argc, argv);

    auto qds(make_shared<DotVector>(3));
    (*qds)[0] = QuantumDot(Eigen::Vector3d(0, 0, 0), 2278.9013,
                           std::pair<double, double>(10, 10),
                           Eigen::Vector3d(5.2917721e-4, 0, 0));
    (*qds)[1] = QuantumDot(Eigen::Vector3d(0, 0, 0.01), 2278.9013,
                           std::pair<double, double>(10, 10),
                           Eigen::Vector3d(5.2917721e-4, 0, 0));
    (*qds)[2] = QuantumDot(Eigen::Vector3d(0, 0.01, 0.01), 2278.9013,
                           std::pair<double, double>(10, 10),
                           Eigen::Vector3d(5.2917721e-4, 0, 0));

    // Set up a dummy Gaussian history
    auto history(make_shared<History::HistoryArray>(boost::extents[3][1000][2]));
    std::fill(history->data(), history->data() + history->num_elements(), Eigen::Vector2cd::Zero());
    for(int i = 0; i < 1000; ++i) {
      const double time = i * config.dt;
      (*history)[0][i][0] =
          Eigen::Vector2cd(0, gaussian((time - config.simulation_time / 2) /
                                       (config.simulation_time / 10)));
    }

    ofstream ipol("input_polarization.dat");
    ipol << scientific << setprecision(12);
    for(int i = 0; i < 1000; ++i) {
      ipol << i * config.dt << " " << polarization((*history)[0][i][0]) << endl;
    }

    // Initalize an interaction
    GreenFunction::RotatingDyadic dy(config.mu0, config.c0, config.omega);
    HistoryInteraction history_interaction(qds, history,
                                           config.interpolation_order, dy);

    ofstream fd("data.dat");
    fd << scientific << setprecision(12);
    for(int i = 0; i < 1000; ++i) {
      history_interaction.evaluate(i);

      fd << i * config.dt << " " << history_interaction.result(0) << " "
         << history_interaction.result(1) << " "
         << history_interaction.result(2) << endl;
    }

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
