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

    vector<QuantumDot> qds(1);
    qds[0] =
        QuantumDot(Eigen::Vector3d(0, 0, 0), 2278.9,
                   std::pair<double, double>(10, 20), Eigen::Vector3d(0, 0, 1));

    cout << setprecision(12) << scientific;
    const double dt = 0.000137856;
    PredictorCorrector::Integrator pc(qds, 7254, dt, 18, 22, 3.15, 1e-12);
    for(int i = 0; i < 7254; ++i) {
      cout << i*dt << " ";
      cout << pc.history[0][i][0][0].real() << " ";
      cout << pc.history[0][i][0][1].real() << " ";
      cout << pc.history[0][i][0][1].imag() << " ";
      cout << endl;

      pc.step();
    }

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
