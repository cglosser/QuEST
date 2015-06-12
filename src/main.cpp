#include <iostream>
#include <iterator>
#include <vector>
#include "input.h"
#include "universe.h"
using namespace std;

int main(int argc, char *argv[]) {
  try {
    auto vm = parse_configs(argc, argv);
    populate_universe(vm);

    cout << "speed of light: " << Universe.c0 << endl;
    cout << "          hbar: " << Universe.hbar << endl;
  } catch(exception &e) {
    cerr << e.what() << endl;
  }
  return 0;
}
