#include <iostream>
#include <vector>

#include "configuration.h"
#include "interaction.h"
#include "lagrange_set.h"
#include "math_utils.h"
#include "quantum_dot.h"
using namespace std;

int main(int argc, char *argv[]) {
  try {
    auto vm = parse_configs(argc, argv);

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
