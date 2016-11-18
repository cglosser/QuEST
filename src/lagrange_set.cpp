#include "lagrange_set.h"

UniformLagrangeSet::UniformLagrangeSet(const int n)
    : order(n), weights(boost::extents[3][n + 1]) {}

UniformLagrangeSet::UniformLagrangeSet(const int n, const double x)
    : UniformLagrangeSet(n)
{
  assert(x >= 0); //Can only interpolate past (i.e. known) values
  calculate_weights(x);
}

void UniformLagrangeSet::calculate_weights(const double x)
{
  for(int basis_id = 0; basis_id <= order; ++basis_id) {
    double d0_product = 1, d1_sum = 0, d2_sum = 0;
    for(int m = 0; m <= order; ++m) {
      if(m == basis_id) continue;
      d0_product *= (x - m)/(basis_id - m);
      d1_sum -= 1/(x - m);
      d2_sum -= std::pow(x - m, -2);
    }

    const size_t idx = basis_id;
    weights[0][idx] = d0_product;
    weights[1][idx] = weights[0][idx]*d1_sum;
    weights[2][idx] = weights[0][idx]*d2_sum +
                      std::pow(weights[1][idx], 2)/weights[0][idx];
  }
}
