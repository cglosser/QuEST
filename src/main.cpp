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

double rhs(const double, const double);
double sol(const double);

int main(int argc, char *argv[]) {
  try {
    auto vm = parse_configs(argc, argv);

    double dt;
    const int WINDOW = 22;
    const vector<double> times(linspace(0, 20, 1024, &dt));

    cout.precision(12);
    cout << scientific;

    PredictorCorrector pc(18, WINDOW, 3.15, 1e-9);

    const double ALPHA = pc.step_factor*dt;

    cout << pc.predictor_coefs << endl << endl;
    cout << pc.corrector_coefs << endl << endl;
    cout << ALPHA << " " << pc.step_factor << " " << dt << endl;

    vector<double>  solution(1024, 0);
    vector<double> dsolution(1024, 0);

    for(size_t i = 0; i < solution.size(); ++i) {
      if(i < WINDOW) {
        solution[i] = sol(times[i]);
        dsolution[i] = rhs(solution[i], times[i]);
      } else {
        for(int w = 0; w < WINDOW; ++w) {
          solution[i] += pc.predictor_coefs(w, 0)*solution[i - WINDOW + w]
                         + pc.predictor_coefs(w, 1)*dsolution[i - WINDOW + w]*ALPHA;
        }
        dsolution[i] = rhs(solution[i], times[i]);

        for(int corr_id = 0; corr_id < 10; ++corr_id) {
          solution[i] = ALPHA*pc.future_coef*dsolution[i];
          for(size_t w = 0; w < WINDOW; ++w) {
            solution[i] += pc.corrector_coefs(w, 0)*solution[i - WINDOW + w]
                           + pc.corrector_coefs(w, 1)*dsolution[i - WINDOW + w]*ALPHA;
          }
          dsolution[i] = rhs(solution[i], times[i]);
        }
      }
    }

    ofstream output("data.dat");
    output << setprecision(12) << scientific;
    for(size_t i = 0; i < solution.size(); ++i) {
      output << times[i] << " " << sol(times[i]) << " " << solution[i] << " " << dsolution[i] << endl;
    }

  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}

double rhs(const double phi, const double x)
{
  return 1/(1 + exp(-(x - 10))) - phi;
}

double sol(const double x)
{
  return (-1 + exp(x) + exp(10)*log(1 + exp(10)) - exp(10)*log(exp(10) + exp(x)))/exp(x);
}
