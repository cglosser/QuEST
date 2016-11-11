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

    cout << setprecision(18) << scientific;

    cout << pc.predictor_coefs << endl << endl;
    cout << pc.corrector_coefs << endl << endl;

    cout << pc.future_coef << endl;


  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
