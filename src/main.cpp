#include <iostream>
#include <iomanip>
#include <iterator>
#include <vector>

#include "bloch.h"
#include "configuration.h"
#include "lagrange_set.h"
using namespace std;

int main(int argc, char *argv[]) {
  try {
    auto vm = parse_configs(argc, argv);

    cout << "      num_particles: " << config.num_particles       << endl;
    cout << "interpolation_order: " << config.interpolation_order << endl;
    cout << "           duration: " << config.simulation_time     << endl;
    cout << "     speed of light: " << config.c0                  << endl;
    cout << "               hbar: " << config.hbar                << endl;

    UniformLagrangeSet uls(-0.2);
    cout << uls.sample_x << endl;

    for(int r = 0; r < 3; ++r) {
      for(int c = 0; c <= config.interpolation_order; ++c) {
        cout << fixed << setprecision(14) << uls.coefficients[r][c] << " ";
      }
      cout << endl;
    }

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
