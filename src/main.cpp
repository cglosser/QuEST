#include <complex>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <vector>

#include "configuration.h"
#include "convolution.h"
#include "interaction.h"
#include "lagrange_set.h"
#include "math_utils.h"
#include "quantum_dot.h"
using namespace std;

typedef std::complex<double> cmplx;

int main(int argc, char *argv[]) {
  try {
    auto vm = parse_configs(argc, argv);

    cout << "      num_particles: " << config.num_particles       << endl;
    cout << "interpolation_order: " << config.interpolation_order << endl;
    cout << "           duration: " << config.simulation_time     << endl;
    cout << "      timestep (dt): " << config.dt                  << endl;
    cout << "     speed of light: " << config.c0                  << endl;
    cout << "               hbar: " << config.hbar                << endl;

    Eigen::Vector3d r1(0,0,0);
    QuantumDot src(Eigen::Vector3d(0,0,0), 0, std::pair<double, double>(0, 0),
                   1, Eigen::Vector3d(0,0,1));

    QuantumDot obs(Eigen::Vector3d(0,0,40.5), 0, std::pair<double, double>(0, 0),
                   1, Eigen::Vector3d(0,0,1));

    cout << "Source: " << src << endl;
    cout << "   Obs: " << obs << endl;

    const cmplx rho00(0, 0);

    ofstream polar("polarization.dat", ios::out),
             efield("efield.dat", ios::out);

    for(int time = 0; time < 45; ++time) {
      cmplx input(gaussian((time - 1024)/128.0), 0);

      src.history.push_back(Eigen::Vector2cd(rho00, input));
      obs.history.push_back(Eigen::Vector2cd(rho00, rho00));

      polar << time << " " << src.polarization(time).transpose() << endl;
      efield << time << " " << r1.transpose() << endl;
    }

    Eigen::Vector3d dr = obs.pos - src.pos;
    Interaction inter(dr);
    cout << " Delay: " << inter.delay.first << " " << inter.delay.second << endl;

    for(int time = 45; time < 2048; ++time) {
      cmplx input(gaussian((time - 1024)/128.0), 0);

      src.history.push_back(Eigen::Vector2cd(rho00, input));
      obs.history.push_back(Eigen::Vector2cd(rho00, rho00));

      auto vecs(inter(src, obs));

      polar << time << " " << src.polarization(time).transpose() << endl;
      efield << time << " " << vecs.second.transpose() << endl;
    }


  } catch(CommandLineException &e) {
    // User most likely queried for help or version info, so we can silently
    // move on
  }

  return 0;
}
