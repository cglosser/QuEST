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

    vector<QuantumDot> dots(2);
    dots[0].pos = Eigen::Vector3d(0, 0, 0);
    dots[0].dipole = Eigen::Vector3d(2, 1, 7);

    dots[1].pos = Eigen::Vector3d(1, 1, 0);
    dots[1].dipole = Eigen::Vector3d(6, 0, -6);

    InteractionTable itab(config.interpolation_order, dots);

    // Build source currents
    for(int i = 0; i < 2048; ++i) {
      Eigen::Vector2cd val(0, gaussian((i - 1024) / 256.0));
      dots[0].history.push_back(val);
    }

    // Evaluate observed "field"
    for(int i = config.interpolation_order + 1; i < 2048; ++i) {
      double val = 0;
      for(int j = 0; j <= config.interpolation_order; ++j) {
        val += 2 * dots[0].history[i - j - itab.floor_delays[0]][1].real() *
               itab.coefficients[0][j];
      }

      cout << setprecision(15) << scientific << i << " "
           << 2 * dots[0].history[i][1].real() << " " << val << endl;
    }

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
