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


    auto in(input_signal());
    auto out(output_signal());

    for(int i = 0; i < 101; ++i) {
      cout << i << "\t" << in[i] << "\t" << out[i] << endl;
    }

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
