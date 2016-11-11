#include <complex>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <vector>

#include "configuration.h"
#include "integrator.h"
#include "interaction.h"
#include "lagrange_set.h"
#include "math_utils.h"
#include "quantum_dot.h"
using namespace std;

typedef std::complex<double> cmplx;

int main(int argc, char *argv[]) {
  try {
    auto vm = parse_configs(argc, argv);

    cout << "      num_particles: " << config.num_particles       << endl;
    cout << "interpolation_order: " << config.interpolation_order << endl;
    cout << "           duration: " << config.simulation_time     << endl;
    cout << "      timestep (dt): " << config.dt                  << endl;
    cout << "     speed of light: " << config.c0                  << endl;
    cout << "               hbar: " << config.hbar                << endl;

    PredictorCorrector pc(18, 22, 3.15, 1e-9);

    cout << pc.corrector_matrix().col(2*pc.n_time) << endl;



  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
