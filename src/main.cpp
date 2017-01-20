#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

#include "configuration.h"
#include "interactions/history_interaction.h"
#include "interactions/rotating_green_function.h"
#include "interactions/history_interaction.h"
#include "math_utils.h"

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
    (*qds)[1] = QuantumDot(Eigen::Vector3d(0, 0, 0.01), 2278.9013,
                           std::pair<double, double>(10, 10),
                           Eigen::Vector3d(5.2917721e-4, 0, 0));
    (*qds)[2] = QuantumDot(Eigen::Vector3d(0, 0.01, 0.01), 2278.9013,
                           std::pair<double, double>(10, 10),
                           Eigen::Vector3d(5.2917721e-4, 0, 0));

    // Set up a dummy Gaussian history
    ofstream initialp("polarization.dat");
    initialp << scientific << setprecision(12);
    auto history(make_shared<History::HistoryArray>(boost::extents[config.num_particles][1000][2]));
    std::fill(history->data(), history->data() + history->num_elements(), Eigen::Vector2cd::Zero());
    for(int i = 0; i < 1000; ++i) {
      const double time = i * config.dt;
      (*history)[0][i][0] =
          Eigen::Vector2cd(0, gaussian((time - config.simulation_time / 2) /
                                       (config.simulation_time / 10)));

      (*history)[1][i][0] = Eigen::Vector2cd(0, sin(M_PI * time / 2));

      initialp << time << " ";
      for(int j = 0; j < config.num_particles; ++j) {
        initialp << (*history)[j][i][0][1] << " ";
      }
      initialp << endl;
    }

    auto gf = std::static_pointer_cast<GreenFunction::Dyadic>(
        std::make_shared<GreenFunction::RotatingDyadic>(config.mu0, config.c0,
                                                        config.omega));

    HistoryInteraction history_interaction(qds, history, gf,
                                           config.interpolation_order);

    ofstream fd("data.dat");
    fd << scientific << setprecision(12);
    for(int i = 0; i < 1000; ++i) {
      history_interaction.evaluate(i);

      fd << i * config.dt << " ";
      for(int j = 0; j < config.num_particles; ++j) {
        fd << history_interaction.result(j) << " ";
      }
      fd << endl;
    }

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
