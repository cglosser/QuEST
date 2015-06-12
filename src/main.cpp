#include <iostream>
#include <iterator>
#include <vector>
#include "input.h"
using namespace std;

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream &os, const vector<T> &v) {
  copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
  return os;
}

int main(int argc, char *argv[]) {
  try {
    auto vm = parse_configs(argc, argv);

    cout << "speed of light: " << vm["constants.c0"].as<double>()   << endl;
    cout << "          hbar: " << vm["constants.hbar"].as<double>() << endl;
  } catch(exception &e) {
    cerr << e.what() << endl;
  }
  return 0;
}
