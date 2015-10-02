#include <iostream>
#include <iterator>
#include <vector>
#include "EMRG_exceptions.h"
#include "input.h"
#include "universe.h"
#include "prolates.h"
using namespace std;

int main(int argc, char *argv[]) {
  try {
    auto vm = parse_configs(argc, argv);

    cout << "    dimensions: " << Universe.dimensions << endl;
    cout << " num_particles: " << Universe.num_particles << endl;
    cout << "      duration: " << Universe.simulation_time << endl;
    cout << "speed of light: " << Universe.c0 << endl;
    cout << "          hbar: " << Universe.hbar << endl;

    Prolate p0(5);

    cout << p0.d0(0) << "\t" << p0.d0(2.2) << "\t" << p0.d0(3) << endl;
    cout << p0.d1(1e-8) << "\t" << p0.d1(2.2) << "\t" << p0.d1(3) << endl;
    cout << p0.d2(1e-4) << "\t" << p0.d2(2.2) << "\t" << p0.d2(3) << endl;
    cout << p0.d3(1e-4) << "\t" << p0.d3(2.2) << "\t" << p0.d3(3) << endl;

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can just bail out
    return 0;
  }
  return 0;
}
