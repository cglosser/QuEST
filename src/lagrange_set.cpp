#include "lagrange_set.h"

Interpolation::UniformLagrangeSet::UniformLagrangeSet(const int order)
    : evaluations(boost::extents[Interpolation::NUM_DERIVATIVES][order + 1]),
      order_(order)
{
}

Interpolation::UniformLagrangeSet::UniformLagrangeSet(const double x,
                                                      const int order,
                                                      const double dt /* = 1 */)
    : UniformLagrangeSet(order)
{
  assert(x > 0);  // Don't extrapolate!
  evaluate_derivative_table_at_x(x, dt);
}

void Interpolation::UniformLagrangeSet::evaluate_derivative_table_at_x(
    const double x, const double dt /* = 1 */)
{
  for(int basis_id = 0; basis_id <= order_; ++basis_id) {
    double d0_product = 1, d1_sum = 0, d2_sum = 0;
    for(int m = 0; m <= order_; ++m) {
      if(m == basis_id) continue;

      double numer = (x - m);

      d0_product *= numer / (basis_id - m);

      if(numer != 0) {
        d1_sum -= 1 / numer;  // Note the minus sign!
        d2_sum -= std::pow(numer, -2);
      }
    }

    evaluations[0][basis_id] = d0_product;
    evaluations[1][basis_id] = (evaluations[0][basis_id] * d1_sum);
    evaluations[2][basis_id] =
        (evaluations[0][basis_id] * d2_sum +
         std::pow(evaluations[1][basis_id], 2) / evaluations[0][basis_id]);
  }

  for(int i = 0; i <= order_; ++i) {
    evaluations[1][i] *= std::pow(dt, -1);
    evaluations[2][i] *= std::pow(dt, -2);
  }
}

Interpolation::DerivFive::DerivFive(const double dt)
    : evaluations(boost::extents[Interpolation::NUM_DERIVATIVES][5 + 1]),
      dt_{dt}
{
}

void Interpolation::DerivFive::evaluate_derivative_table_at_x(const double x)
{
  for(int basis_id = 0; basis_id <= 5; ++basis_id) {
    evaluations[0][basis_id] = d0(basis_id, x);
    evaluations[1][basis_id] = d1(basis_id, x) / std::pow(dt_, 1);
    evaluations[2][basis_id] = d2(basis_id, x) / std::pow(dt_, 2);
  }
}

// clang-format off
double Interpolation::DerivFive::d0(const int i, const double x) const
{
  switch(i) {
    case 0: return 1+x*(137.0/60+x*(15.0/8+x*(17.0/24+x*(1.0/8+x/120))));
    case 1: return x*(-5+x*(-77.0/12+x*(-71.0/24+x*(-7.0/12-x/24))));
    case 2: return x*(5+x*(107.0/12+x*(59.0/12+x*(13.0/12+x/12))));
    case 3: return x*(-10.0/3+x*(-13.0/2+x*(-49.0/12+x*(-1-x/12))));
    case 4: return x*(5.0/4+x*(61.0/24+x*(41.0/24+x*(11.0/24+x/24))));
    case 5: return x*(-1.0/5+x*(-5.0/12+x*(-7.0/24+x*(-1.0/12-x/120))));
    default: return 0;
  }
}

double Interpolation::DerivFive::d1(const int i, const double x) const
{
  switch(i) {
    case 0: return 137.0/60+x*(15.0/4+x*(17.0/8+x*(1.0/2+x/24)));
    case 1: return -5+x*(-77.0/6+x*(-71.0/8+x*(-7.0/3-(5*x)/24)));
    case 2: return 5+x*(107.0/6+x*(59.0/4+x*(13.0/3+(5*x)/12)));
    case 3: return -10.0/3+x*(-13+x*(-49.0/4+x*(-4-(5*x)/12)));
    case 4: return 5.0/4+x*(61.0/12+x*(41.0/8+x*(11.0/6+(5*x)/24)));
    case 5: return -1.0/5+x*(-5.0/6+x*(-7.0/8+x*(-1.0/3-x/24)));
    default: return 0;
  }
}

double Interpolation::DerivFive::d2(const int i, const double x) const
{
  switch(i) {
    case 0: return 15.0/4+x*(17.0/4+x*(3.0/2+x/6));
    case 1: return -77.0/6+x*(-71.0/4+x*(-7-(5*x)/6));
    case 2: return 107.0/6+x*(59.0/2+x*(13+(5*x)/3));
    case 3: return -13+x*(-49.0/2+x*(-12-(5*x)/3));
    case 4: return 61.0/12+x*(41.0/4+x*(11.0/2+(5*x)/6));
    case 5: return -5.0/6+x*(-7.0/4+x*(-1-x/6));
    default: return 0;
  }
}
// clang-format on
