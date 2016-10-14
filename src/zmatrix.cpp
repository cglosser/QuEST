#include "zmatrix.h"

std::vector<double> lagrange_coefficients(
    const int interp_order,
    const int deriv_order,
    const double x)
{
  assert(x <= 0); //Can only interpolate past (i.e. known) values
  std::vector<double> coefficients(interp_order + 1, 1);

  for(int basis_id = -interp_order; basis_id <= 0; ++basis_id) {
    for(int m = -interp_order; m <= 0; ++m) {
      if(m == basis_id) continue;
      coefficients.at(basis_id + interp_order) *= (m - x)/(m - basis_id);
    }
  }

  return coefficients;
}

zmatrix::zmatrix(const std::vector<Eigen::Vector3d> &pts, 
    const int history_length) :
  weights(boost::extents[pts.size()][pts.size()][history_length])
{

}
