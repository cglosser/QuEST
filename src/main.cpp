#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

#include "bloch_system.h"
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

    vector<QuantumDot> qds(2);
    qds[0] = QuantumDot(Eigen::Vector3d(), 0, std::pair<double, double>(10, 20),
                        Eigen::Vector3d(0, 0, 1));

    qds[1] =
        QuantumDot(Eigen::Vector3d(0.005, 0.005, 0.005), 0,
                   std::pair<double, double>(10, 20), Eigen::Vector3d(0, 0, 1));

    //PredictorCorrector::Integrator rpc(30, 18, 22, 3.15, 1e-12);

    //cout << setprecision(12) << scientific;
    //cout << rpc.weights_.ps << endl << endl;
    //cout << rpc.weights_.cs << endl << endl;


  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
