#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

#include "configuration.h"
#include "interaction_table.h"
#include "lagrange_set.h"
#include "math_utils.h"
#include "quantum_dot.h"

using namespace std;

int main(int argc, char *argv[])
{
  try {
    auto vm = parse_configs(argc, argv);

    DotTable qds(2);

    qds[0] =
        QuantumDot(Eigen::Vector3d(0, 0, 0), 2278.9,
                   std::pair<double, double>(10, 20), Eigen::Vector3d(0, 0, 1));
    qds[1] =
        QuantumDot(Eigen::Vector3d(0.01, 0.01, 0.01), 2278.9,
                   std::pair<double, double>(10, 20), Eigen::Vector3d(0, 0, 1));

    InteractionTable itable(config.interpolation_order, std::move(qds));

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
