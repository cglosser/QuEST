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

    shared_ptr<vector<QuantumDot>> qds(new vector<QuantumDot>(2));
    (*qds)[0] = QuantumDot(Eigen::Vector3d(0, 0, 0), 100,
                           std::pair<double, double>(10, 20),
                           Eigen::Vector3d(1.00000, 0, 0));
    (*qds)[1] = QuantumDot(Eigen::Vector3d(0.1, 0.2, 0.3), 100,
                           std::pair<double, double>(10, 20),
                           Eigen::Vector3d(3.14159, 0, 0));

    std::shared_ptr<const Pulse> p(new Pulse(
        1, 0, 0.1, 1, Eigen::Vector3d(0, 0, 1), Eigen::Vector3d(1, 0, 0)));

    InteractionTable it(3, qds, p);

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
