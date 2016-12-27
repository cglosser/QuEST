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
    auto vm = parse_configs(argc, argv);

    shared_ptr<DotVector> dots(new DotVector(import_dots(config.qd_path)));
    shared_ptr<Pulse> pulse(new Pulse(1558.92260227, 5, 2278.9013, 2278.9013,
                                      Eigen::Vector3d(0, 0, 1),
                                      Eigen::Vector3d(1, 0, 0)));
    InteractionTable interaction_table(config.interpolation_order, dots, pulse);
    auto superoperators(rhs_functions(*dots));

    const int num_steps =
        static_cast<int>(ceil(config.simulation_time / config.dt));
    PredictorCorrector::Integrator solver(2, num_steps, config.dt, 18, 22, 3.15,
                                          std::move(superoperators),
                                          interaction_table);

    ofstream bloch("bloch.dat");
    bloch << scientific << setprecision(14);

    for(int i = 1; i < num_steps; ++i) {
      solver.solve(i);
      bloch << i*config.dt << " " << solver.history[0][i][0].transpose() << " " << solver.history[1][i][0].transpose() << endl;
    }

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
