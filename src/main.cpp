#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

#include "configuration.h"
#include "interactions/green_function.h"
#include "interactions/history_interaction.h"
#include "quantum_dot.h"

using namespace std;

int main(int argc, char *argv[])
{
  try {
    cout << setprecision(12) << scientific;
    auto vm = parse_configs(argc, argv);


    auto qds(make_shared<DotVector>(2));
    (*qds)[0] =
        QuantumDot(Eigen::Vector3d(0, 0, 0), 2278.9013,
                   std::pair<double, double>(10, 10), Eigen::Vector3d(5.2917721e-4, 0, 0));
    (*qds)[1] =
        QuantumDot(Eigen::Vector3d(0, 0, 0.01), 2278.9013,
                   std::pair<double, double>(10, 10), Eigen::Vector3d(5.2917721e-4, 0, 0));

    auto hist(make_shared<History::HistoryArray>());

    auto ffd(make_shared<GreenFunction::Dyadic>(config.mu0, config.c0));

    HistoryInteraction(qds, hist, ffd, config.interpolation_order);

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
