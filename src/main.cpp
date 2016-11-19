#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

#include "bloch_system.h"
#include "configuration.h"
#include "interaction_table.h"
#include "lagrange_set.h"
#include "math_utils.h"
#include "quantum_dot.h"

using namespace std;

int main(int argc, char *argv[])
{
  try {
    auto vm = parse_configs(argc, argv);

    BlochSystem bs(import_dots("dots.dat"));

    cout << bs.dots[0].history.capacity() << endl;


  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
