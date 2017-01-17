#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

#include "configuration.h"
/*
 *#include "integrator.h"
 *#include "interactions/history_interaction.h"
 *#include "lagrange_set.h"
 *#include "math_utils.h"
 *#include "pulse.h"
 *#include "quantum_dot.h"
 */

#include "interactions/green_function.h"

using namespace std;

int main(int argc, char *argv[])
{
  try {
    cout << setprecision(12);
    auto vm = parse_configs(argc, argv);

    GreenFunction::Dyadic ffd(1, 1);
    Interpolation::UniformLagrangeSet uls(0.5, 3);

    auto dyads(ffd.calculate_coefficients(Eigen::Vector3d(0.1, 0.2, 0.3), uls));

    for(const auto &dyad : dyads) { cout << dyad.real() << endl << endl; }


  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
