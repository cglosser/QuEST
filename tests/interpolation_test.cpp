#include <cmath>
#include <iostream>
#include <vector>

#include <boost/math/special_functions/bessel.hpp>

#include "prolates.h"

using namespace std;

#define WIDTH 4
#define NSAMP 16

double hat(double x) {
  if (x > 0) {
    return 1 - x;
  } else if (x <= 0) {
    return x + 1;
  } else {
    return 0;
  }
}

int main()
{

  const double t_max = boost::math::cyl_bessel_j_zero(0.0, 5);
  const double dt = 1;//t_max/(NSAMP - 1);
  const size_t nsteps = t_max/dt;

  double jnaught = boost::math::cyl_bessel_j(0.0, 0.0);
  ProlateTimeExpansion pte(Prolate(3), nsteps + 5, Eigen::Vector3d(jnaught, jnaught, jnaught));

  for(int i = 1; i < nsteps; ++i) {
    double time = i*dt;
    double val = boost::math::cyl_bessel_j(0.0, time);

    pte.step(Eigen::Vector3d(val, val, val));
  }

  for(double time = 0; time < t_max; time += 0.001) {
    cout << time << "\t" << pte.at(time/dt)[0] << endl;
  }

}
