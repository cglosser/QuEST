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
#include "quantum_dot.h"

using namespace std;

int main(int argc, char *argv[])
{
  try {
    cout << setprecision(12) << scientific;
    auto vm = parse_configs(argc, argv);

    DotTable qds(2);

    const double dt = config.dt;
    const int nsteps = 1/dt;
    using Hist = PredictorCorrector::Integrator::HistoryArray;

    qds[0] =
        QuantumDot(Eigen::Vector3d(0, 0, 0), 2278.9,
                   std::pair<double, double>(10, 20), Eigen::Vector3d(0, 0, 1));
    qds[1] =
        QuantumDot(Eigen::Vector3d(0.01, 0.01, 0.01), 2278.9,
                   std::pair<double, double>(10, 20), Eigen::Vector3d(0, 0, 1));

    InteractionTable itable(config.interpolation_order, qds);
    Hist ha(boost::extents[qds.size()][Hist::extent_range(-22, nsteps)][2]);

    for(int i = -22; i < nsteps; ++i) {
      ha[0][i][0] = matrix_elements(0, gaussian((i*dt - 0.5)/0.1));
    }

    for(int i = 0; i < nsteps; ++i) {
      itable.convolve_currents(ha, i);
      cout << i << " " << i*dt << " " << polarization(ha[0][i][0]) << " " << itable.convolution.at(1) << endl;
    }
  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
