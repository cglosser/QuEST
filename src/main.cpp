#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

#include "configuration.h"
#include "integrator.h"
#include "interaction_table.h"
#include "lagrange_set.h"
#include "math_utils.h"
#include "pulse.h"
#include "quantum_dot.h"

using namespace std;

int main(int argc, char *argv[])
{
  try {
    cout << setprecision(12) << scientific;
    auto vm = parse_configs(argc, argv);

    Pulse pulse(1, 2, 0.2, 1, Eigen::Vector3d(0,0,1), Eigen::Vector3d(2,3,4));
    Pulse p;

    cin >> p;
    cout << p;

    //for(int i = -1000; i <=1000; ++i) {
      //const double x = i/250.;

      //cout << x << " " << pulse(Eigen::Vector3d(0,0,x), 0).transpose() << endl;
    //}

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
