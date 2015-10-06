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

  vector<double> weights;
  for(int i = -WIDTH; i < NSAMP + WIDTH; ++i) {
    double val = boost::math::cyl_bessel_j(0.0, i*dt);
    weights.push_back(val);
  }

  Prolate p0(WIDTH);

  double step = dt/8;
  for(double x = 0; x < t_max; x += step) {
    int base_idx = floor(x/dt + 0.5);

    double val = 0;
    for(int w = -WIDTH; w < WIDTH; ++w) {
      val += weights[base_idx + w + WIDTH]*p0.d0(x/dt - base_idx - w);
    }

    cout << x << "\t" << "\t" << boost::math::cyl_bessel_j(0.0, x) << "\t" << val << endl;
  }

}
