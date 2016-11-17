#include <complex>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <vector>

#include "configuration.h"
#include "interaction.h"
#include "lagrange_set.h"
#include "math_utils.h"
#include "quantum_dot.h"
using namespace std;

typedef std::complex<double> cmplx;

int main(int argc, char *argv[]) {
  try {
    auto vm = parse_configs(argc, argv);

    QuantumDot dot1(Eigen::Vector3d(0,0,0), 1000, std::pair<double, double>(10, 10), Eigen::Vector3d(0,0,1));
    QuantumDot dot2(Eigen::Vector3d(8,8,3.1415926535), 1000, std::pair<double, double>(10, 10), Eigen::Vector3d(0,0,1));

    Interaction inter(dot1, dot2);

    cout << inter.farfield_dyadic(dot2.pos - dot1.pos, dot1.dipole, dot2.dipole) << endl;
    cout << inter.midfield_dyadic(dot2.pos - dot1.pos, dot1.dipole, dot2.dipole) << endl;
    cout << inter.nearfield_dyadic(dot2.pos - dot1.pos, dot1.dipole, dot2.dipole) << endl;

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
