#include <iostream>
#include <iterator>
#include <vector>

#include "bloch.h"
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

    ProlateTimeExpansion pte(Prolate(2), 10, Eigen::Vector3d(10,10,10));

    for(int i = 0; i < 10; ++i) {
      pte.step(Eigen::Vector3d(10,10,10));
    }

    for(int i = 0; i <= 4; ++i) {
      cout << pte.at(i + 0.2).transpose() << endl;
    }


  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can just bail out
    return 0;
  }
  return 0;
}
