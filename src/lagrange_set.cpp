#include "lagrange_set.h"

UniformLagrangeSet::UniformLagrangeSet(const double x, const int ord) :
  sample_x(x), order(ord), coefficients(boost::extents[3][ord + 1])
{
  assert(x <= 0); //Can only interpolate past (i.e. known) values
  set_deriv_0();
  set_deriv_1();
  set_deriv_2();
}

void UniformLagrangeSet::set_deriv_0()
{
  for(int basis_id = -order; basis_id <= 0; ++basis_id) {
    double rhs_product = 1;
    for(int m = -order; m <= 0; ++m) {
      if(m == basis_id) continue;
      rhs_product *= (m - sample_x)/(m - basis_id);
    }

    const size_t idx = basis_id + order;
    coefficients[0][idx] = rhs_product;
  }
}

void UniformLagrangeSet::set_deriv_1()
{
  for(int basis_id = -order; basis_id <= 0; ++basis_id) {
    double rhs_sum = 0;
    for(int m = -order; m <= 0; ++m) {
      if(m == basis_id) continue;
      rhs_sum += 1/(sample_x - m);
    }

    const size_t idx = basis_id + order;
    coefficients[1][idx] = coefficients[0][idx]*rhs_sum;
  }
}

void UniformLagrangeSet::set_deriv_2()
{
  for(int basis_id = -order; basis_id <=0; ++basis_id) {
    double rhs_sum = 0;
    for(int m = -order; m <= 0; ++m) {
      if(m == basis_id) continue;
      rhs_sum -= std::pow(sample_x - m, -2);
    }

    const size_t idx = basis_id + order;
    coefficients[2][idx] =
      coefficients[0][idx]*rhs_sum
      + std::pow(coefficients[1][idx], 2)/coefficients[0][idx];
  }
}