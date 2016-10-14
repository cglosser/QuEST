#include <iostream>
#include <iterator>
#include <vector>

#include "bloch.h"
#include "input.h"
#include "universe.h"
#include "prolates.h"
#include "zmatrix.h"
using namespace std;

int main(int argc, char *argv[]) {
  try {
    auto vm = parse_configs(argc, argv);

    cout << "    dimensions: " << Universe.dimensions << endl;
    cout << " num_particles: " << Universe.num_particles << endl;
    cout << "      duration: " << Universe.simulation_time << endl;
    cout << "speed of light: " << Universe.c0 << endl;
    cout << "          hbar: " << Universe.hbar << endl;

    //vector<Eigen::Vector3d> pts(42);

    //zmatrix myMat(pts, 816);

    for(auto it : lagrange_coefficients(3, -0.2)) {
      cout << it << " ";
    }
    cout << endl;

    for(auto it : deriv_lagrange_coefficients(3, -0.2)) {
      cout << it << " ";
    }
    cout << endl;



  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can just bail out
    return 0;
  }
  return 0;
}
