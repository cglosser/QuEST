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
    qds[0].pos = Eigen::Vector3d(0, 0, 0);
    qds[0].dipole = Eigen::Vector3d(0, 0, 1);

    qds[1].pos = Eigen::Vector3d(16, 8, 4);
    qds[1].dipole = Eigen::Vector3d(0, 0, 1);

    PredictorCorrector rpc(18, 22, 3.15, 1e-12);

    BlochSystem sys(rpc, qds, config.interpolation_order, 2048);

    for(int i = -22; i < 2048; ++i) {
      const double g = skew_gaussian((i - 1024)/256., 5);

      sys.history[0][i][0][1] = g;
      sys.history[1][i][0][1] = 0;
    }

    for(int i = 0; i < 2048; ++i) {
      sys.convolve_currents();
      sys.now++;

      cout << setprecision(12);
      cout << i << " " << sys.rabi_freqs[0] << " " << sys.rabi_freqs[1] << endl;
    }

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
