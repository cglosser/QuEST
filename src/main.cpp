#include <iostream>
#include <iterator>
#include <vector>
#include "EMRG_exceptions.h"
#include "input.h"
#include "universe.h"
using namespace std;

int main(int argc, char *argv[]) {
  try {
    auto vm = parse_configs(argc, argv);

    cout << "    dimensions: " << Universe.dimensions << endl;
    cout << " num_particles: " << Universe.num_particles << endl;
    cout << "      duration: " << Universe.simulation_time << endl;
    cout << "speed of light: " << Universe.c0 << endl;
    cout << "          hbar: " << Universe.hbar << endl;
  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can just bail out
    return 0; 
  }
  return 0;
}
