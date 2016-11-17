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

    QuantumDot dot1(Eigen::Vector3d(0,0,0), 1000, std::pair<double, double>(10, 10), 1, Eigen::Vector3d(1,2,3));
    QuantumDot dot2(Eigen::Vector3d(3,4,5), 1000, std::pair<double, double>(10, 10), 1, Eigen::Vector3d(1,2,3));

    Interaction inter(dot1, dot2);

    cout << inter.delay.first << " " << inter.delay.second << endl;

    cout << inter.rhat_dyadic(dot1.pos - dot2.pos) << endl;


  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
