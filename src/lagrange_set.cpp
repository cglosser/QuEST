#include "lagrange_set.h"


Interpolation::UniformLagrangeSet::UniformLagrangeSet(const int order) :
  weights(boost::extents[Interpolation::NUM_DERIVATIVES][order + 1]),
  order_(order) { }

Interpolation::UniformLagrangeSet::UniformLagrangeSet(const double x,
                                                      const int order,
                                                      const double dt /* = 1 */)
    : UniformLagrangeSet(order)
{
  assert(x > 0);  // Don't extrapolate!
  calculate_weights(x, dt);
}

void Interpolation::UniformLagrangeSet::calculate_weights(
    const double x, const double dt /* = 1 */)
{
  for(int basis_id = 0; basis_id <= order_; ++basis_id) {
    double d0_product = 1, d1_sum = 0, d2_sum = 0;
    for(int m = 0; m <= order_; ++m) {
      if(m == basis_id) continue;
      d0_product *= (x - m) / (basis_id - m);
      d1_sum -= 1 / (x - m); // Note the minus sign!
      d2_sum -= std::pow(x - m, -2);
    }

    weights[0][basis_id] = d0_product;
    weights[1][basis_id] = (weights[0][basis_id] * d1_sum);
    weights[2][basis_id] =
        (weights[0][basis_id] * d2_sum +
         std::pow(weights[1][basis_id], 2) / weights[0][basis_id]);
  }

  for(int i = 0; i <= order_; ++i) {
    weights[1][i] *= std::pow(dt, -1);
    weights[2][i] *= std::pow(dt, -2);
  }

}
