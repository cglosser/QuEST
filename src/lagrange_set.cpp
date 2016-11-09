#include "lagrange_set.h"

UniformLagrangeSet::UniformLagrangeSet()
    : sample_x(0), weights(boost::extents[3][config.interpolation_order + 1])
{
  for(int r = 0; r < 3; ++r) {
    for(int c = 0; c <= config.interpolation_order; ++c) {
      weights[r][c] = 0;
    }
  }
}

UniformLagrangeSet::UniformLagrangeSet(const double x)
    : sample_x(x), weights(boost::extents[3][config.interpolation_order + 1])
{
  assert(x >= 0); //Can only interpolate past (i.e. known) values
  calculate_weights();
}

void UniformLagrangeSet::calculate_weights()
{
  for(int basis_id = 0; basis_id <= config.interpolation_order; ++basis_id) {
    double d0_product = 1, d1_sum = 0, d2_sum = 0;
    for(int m = 0; m <= config.interpolation_order; ++m) {
      if(m == basis_id) continue;
      d0_product *= (sample_x - m)/(basis_id - m);
      d1_sum -= 1/(sample_x - m); /* I have _no_ idea why this one needs -= (instead of +=) to produce the appropriate
                                   * derivative, but until I can figure it out I'm leaving this comment as a warning to
                                   * those who look at this later.
                                   */
      d2_sum -= std::pow(sample_x - m, -2);
    }

    const size_t idx = basis_id;
    weights[0][idx] = d0_product;
    weights[1][idx] = weights[0][idx]*d1_sum;
    weights[2][idx] =
      weights[0][idx]*d2_sum + std::pow(weights[1][idx], 2)/weights[0][idx];
  }
}
