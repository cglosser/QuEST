#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#include "configuration.h"
#include "integrator.h"
#include "interaction.h"
#include "lagrange_set.h"
#include "math_utils.h"
#include "quantum_dot.h"
using namespace std;

int main(int argc, char *argv[]) {
  try {
    auto vm = parse_configs(argc, argv);

    PredictorCorrector pc(18, 22, 3.15, 1e-9);

    auto x(pc.compute_coefficients(pc.predictor_matrix()));

    ofstream mat("mat.dat"), vec("vec.dat");

    mat << pc.predictor_matrix();
    vec << pc.rhs_vector();

    cout << setprecision(18) << scientific;
    cout << pc.compute_coefficients(pc.predictor_matrix()) << endl;


  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
