#include <Eigen/Dense>
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <vector>

#include "configuration.h"
#include "integrator.h"
#include "interactions/history_interaction.h"
#include "lagrange_set.h"
#include "math_utils.h"
#include "pulse.h"
#include "quantum_dot.h"

using namespace std;

void test_history_interaction(int argc, char *argv[])
{
  cout << setprecision(15) << scientific;
  auto vm = parse_configs(argc, argv);
  auto dots = make_shared<DotVector>(import_dots(config.qd_path));

  const int num_steps = 100000;
  auto history =
      make_shared<History::HistoryArray>(boost::extents[2][num_steps][2]);

  cout << "Distance (dimensionless): "
       << ((*dots)[1].position() - (*dots)[0].position()).norm() /
              (config.c0 * config.dt) << endl;

  for(int time_idx = 0; time_idx < num_steps; time_idx++) {
    const double time = time_idx * config.dt;
    const double signal = gaussian((time - 5)/1);
    (*history)[0][time_idx][0] = Eigen::Vector2cd(0, signal);
  }

  HistoryInteraction hint(dots, history, config.interpolation_order);

  for(int time_idx = 10; time_idx < num_steps; time_idx++) {
    const double time = time_idx * config.dt;
    hint.evaluate(time_idx);
    const double val = hint.result(1);
    (*history)[1][time_idx][1] = Eigen::Vector2cd(val, 0);

    cout << time << " " << (*history)[0][time_idx][0][1].real() << " "
         << (*history)[1][time_idx][1][0].real() << endl;
  }
}

int main(int argc, char *argv[])
{
  try {
    auto vm = parse_configs(argc, argv);

    test_history_interaction(argc, argv);
  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
