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

    PredictorCorrector rpc(18, 22, 3.15, 1e-12);
    vector<QuantumDot> qds(1);
    qds[0] = QuantumDot(Eigen::Vector3d(0,0,0), 0.05, std::pair<double, double>(10, 20), Eigen::Vector3d(0, 0, 1));

    BlochSystem sys(rpc, qds, 2048);

    for(int i = 0; i < 2048; ++i) sys.step();

    for(int i = -22; i < 2048; ++i) {
      cout << setw(5) << i << " ";
      cout << setprecision(12) << scientific << sys.history[0][i][0][0].real() << endl;
    }


  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
