#include <complex>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <vector>

#include "configuration.h"
#include "interaction_table.h"
#include "lagrange_set.h"
#include "math_utils.h"
#include "quantum_dot.h"
using namespace std;

typedef std::complex<double> cmplx;

int main(int argc, char *argv[]) {
  try {
    auto vm = parse_configs(argc, argv);

    vector<QuantumDot> dots(2);
    dots[0].pos = Eigen::Vector3d(0, 0, 0);
    dots[0].dipole = Eigen::Vector3d(2, 1, 7);

    dots[1].pos = Eigen::Vector3d(1, 1, 0);
    dots[1].dipole = Eigen::Vector3d(6, 0, -6);

    InteractionTable itab(config.interpolation_order, dots);

    for(int r = 0; r < itab.num_interactions; ++r) {
      cout << r << " | ";
      for(int c = 0; c <= config.interpolation_order; ++c) {
        cout << setprecision(15) << scientific << itab.coefficients[r][c] << " ";
      }
      cout << endl;
    }

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
