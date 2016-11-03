#include "lagrange_set.h"

UniformLagrangeSet::UniformLagrangeSet(const double x) :
  sample_x(x), weights(boost::extents[3][config.interpolation_order + 1])
{
  assert(x <= 0); //Can only interpolate past (i.e. known) values
  set_deriv_0();
  set_deriv_1();
  set_deriv_2();
}

void UniformLagrangeSet::set_deriv_0()
{
  for(int basis_id = -config.interpolation_order; basis_id <= 0; ++basis_id) {
    double rhs_product = 1;
    for(int m = -config.interpolation_order; m <= 0; ++m) {
      if(m == basis_id) continue;
      rhs_product *= (m - sample_x)/(m - basis_id);
    }

    const size_t idx = basis_id + config.interpolation_order;
    weights[0][idx] = rhs_product;
  }
}

void UniformLagrangeSet::set_deriv_1()
{
  for(int basis_id = -config.interpolation_order; basis_id <= 0; ++basis_id) {
    double rhs_sum = 0;
    for(int m = -config.interpolation_order; m <= 0; ++m) {
      if(m == basis_id) continue;
      rhs_sum += 1/(sample_x - m);
    }

    const size_t idx = basis_id + config.interpolation_order;
    weights[1][idx] = weights[0][idx]*rhs_sum;
  }
}

void UniformLagrangeSet::set_deriv_2()
{
  for(int basis_id = -config.interpolation_order; basis_id <=0; ++basis_id) {
    double rhs_sum = 0;
    for(int m = -config.interpolation_order; m <= 0; ++m) {
      if(m == basis_id) continue;
      rhs_sum -= std::pow(sample_x - m, -2);
    }

    const size_t idx = basis_id + config.interpolation_order;
    weights[2][idx] =
      weights[0][idx]*rhs_sum + std::pow(weights[1][idx], 2)/weights[0][idx];
  }
}
