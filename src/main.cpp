#include <iostream>
#include <iomanip>
#include <iterator>
#include <vector>

#include "bloch.h"
#include "configuration.h"
#include "convolution.h"
#include "interaction.h"
#include "lagrange_set.h"
#include "math_utils.h"
using namespace std;

int main(int argc, char *argv[]) {
  try {
    auto vm = parse_configs(argc, argv);

    cout << "      num_particles: " << config.num_particles       << endl;
    cout << "interpolation_order: " << config.interpolation_order << endl;
    cout << "           duration: " << config.simulation_time     << endl;
    cout << "      timestep (dt): " << config.dt                  << endl;
    cout << "     speed of light: " << config.c0                  << endl;
    cout << "               hbar: " << config.hbar                << endl;

    Eigen::Vector3d r1(0,0,0), r2(300,400,0);

    Interaction inter(r2 - r1);

    cout << (r2 - r1).norm()/(config.c0*config.dt) << endl;
    cout << inter.delay.first << " " << inter.delay.second << endl;

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
