#include "zmatrix.h"

constexpr int NUM_DERIVATIVES = 3;

std::vector<double> deriv_0_lagrange_coefficients(
    const int interp_order,
    const double x)
{
  assert(x <= 0); //Can only interpolate past (i.e. known) values
  std::vector<double> result(interp_order + 1);

  for(int basis_id = -interp_order; basis_id <= 0; ++basis_id) {
    double rhs_product = 1;
    for(int m = -interp_order; m <= 0; ++m) {
      if(m == basis_id) continue;
      rhs_product *= (m - x)/(m - basis_id);
    }

    const size_t idx = basis_id + interp_order;
    result.at(idx) = rhs_product;
  }

  return result;
}

std::vector<double> deriv_1_lagrange_coefficients(
    const int interp_order,
    const double x)
{
  assert(x <= 0); //Can only interpolate past (i.e. known) values
  std::vector<double> d0(deriv_0_lagrange_coefficients(interp_order, x));
  std::vector<double> result(interp_order + 1);

  for(int basis_id = -interp_order; basis_id <= 0; ++basis_id) {
    double rhs_sum = 0;
    for(int m = -interp_order; m <= 0; ++m) {
      if(m == basis_id) continue;
      rhs_sum += 1/(x - m);
    }

    const size_t idx = basis_id + interp_order;
    result.at(idx) = d0.at(idx)*rhs_sum;
  }

  return result;
}

std::vector<double> deriv_2_lagrange_coefficients(
    const int interp_order,
    const double x)
{
  assert(x <= 0); //Can only interpolate past (i.e. known) values
  std::vector<double> d0(deriv_0_lagrange_coefficients(interp_order, x));
  std::vector<double> d1(deriv_1_lagrange_coefficients(interp_order, x));
  std::vector<double> result(interp_order + 1);

  for(int basis_id = -interp_order; basis_id <=0; ++basis_id) {
    double rhs_sum = 0;
    for(int m = -interp_order; m <= 0; ++m) {
      if(m == basis_id) continue;
      rhs_sum -= std::pow(x - m, -2);
    }

    const size_t idx = basis_id + interp_order;
    result.at(idx) = d0.at(idx)*rhs_sum + std::pow(d1.at(idx), 2)/d0.at(idx);
  }

  return result;
}

zmatrix::zmatrix(const std::vector<Eigen::Vector3d> &pts,
    const int history_length) :
  weights(boost::extents[pts.size()][pts.size()][history_length][NUM_DERIVATIVES])
{}
