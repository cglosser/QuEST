#include <iostream>
#include <Eigen/Dense>
#include <boost/multi_array.hpp>


class zmatrix
{
 public:
  typedef boost::multi_array<double, 4> array;
  typedef array::index index;

  array weights;

  zmatrix(const std::vector<Eigen::Vector3d> &, const int);

};

std::vector<double> deriv_0_lagrange_coefficients(int, double);
std::vector<double> deriv_1_lagrange_coefficients(int, double);
