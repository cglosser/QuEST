#include <iostream>
#include <iomanip>
#include <iterator>
#include <vector>

#include "bloch.h"
#include "configuration.h"
#include "convolution.h"
#include "lagrange_set.h"
#include "math_utils.h"
using namespace std;

int main(int argc, char *argv[]) {
  try {
    auto vm = parse_configs(argc, argv);

    cout << "      num_particles: " << config.num_particles       << endl;
    cout << "interpolation_order: " << config.interpolation_order << endl;
    cout << "           duration: " << config.simulation_time     << endl;
    cout << "     speed of light: " << config.c0                  << endl;
    cout << "               hbar: " << config.hbar                << endl;

    QuantumDot qd;
    int ipart = 34;
    UniformLagrangeSet delay(0.6);

    for(int i = 0; i < 500; ++i) {
      double val = exp(-std::pow((i - 250)/25.0, 2)/2);
      qd.history.push_back(val);

      if(i > ipart) {
        cout << i << " ";
        cout << qd.history[i] << " ";
        cout << qd.interpolate(delay, ipart) << endl;
      }
    }

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
