#include <complex>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <vector>

#include "configuration.h"
#include "interaction_table.h"
#include "interaction.h"
#include "lagrange_set.h"
#include "math_utils.h"
#include "quantum_dot.h"
using namespace std;

typedef std::complex<double> cmplx;

int main(int argc, char *argv[]) {
  try {
    auto vm = parse_configs(argc, argv);

    vector<QuantumDot> dots(6);
    for(int i = 0; i < 6; ++i) {
      dots[i].pos = Eigen::Vector3d(i, std::sin(i/6.0), 0);
      dots[i].dipole = Eigen::Vector3d(0, 0, 1);
    }

    InteractionTable itab(dots);

    cout << "Interaction table:" << endl;
    for(int i = 0; i < 15; ++i) {
      cout << i << " | ";
      for(int r = 0; r <= config.interpolation_order; ++r) {
        cout << itab.coefficients[i][r] << " ";
      }
      cout << endl;
    }

    cout << endl << endl;

    cout << "Interaction objects:" << endl;
    int idx = 0;
    for(int d1 = 0; d1 < 5; ++d1) {
      for(int d2 = d1 + 1; d2 < 6; ++d2) {
        Interaction inter(dots[d2], dots[d1]);
        cout << idx++ << " | ";
        for(auto const &x : inter.coefs) cout << x << " ";
        cout << endl;
      }
    }


  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
