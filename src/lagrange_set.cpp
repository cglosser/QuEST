#include "lagrange_set.h"

constexpr int NUM_DERIVATIVES = 3;

UniformLagrangeSet::UniformLagrangeSet(const int n, const double timestep = 1)
    : order(n), dt(timestep), weights(boost::extents[NUM_DERIVATIVES][n + 1])
{
}

UniformLagrangeSet::UniformLagrangeSet(const double x, const int n,
                                       const double timestep = 1)
    : UniformLagrangeSet(n, timestep)
{
  assert(x > 0);  // Don't extrapolate!
  calculate_weights(x);
}

// This calculation assumes the x coordinates of each knot correspond to
// *positive* integers -- because the resulting weights generally get used in
// the calculation of retardation effects, odd derivatives need an extra minus
// sign!
void UniformLagrangeSet::calculate_weights(const double x)
{
  for(int basis_id = 0; basis_id <= order; ++basis_id) {
    double d0_product = 1, d1_sum = 0, d2_sum = 0;
    for(int m = 0; m <= order; ++m) {
      if(m == basis_id) continue;
      d0_product *= (x - m) / (basis_id - m);
      d1_sum -= 1 / (x - m);
      d2_sum -= std::pow(x - m, -2);
    }

    weights[0][basis_id] = d0_product;
    weights[1][basis_id] = (weights[0][basis_id] * d1_sum) / dt;
    weights[2][basis_id] =
        (weights[0][basis_id] * d2_sum +
         std::pow(weights[1][basis_id], 2) / weights[0][basis_id]) /
        std::pow(dt, 2);
  }
}
