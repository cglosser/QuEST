#include <iostream>
#include <Eigen/Dense>
#include <boost/multi_array.hpp>


class zmatrix
{
 public:
  typedef boost::multi_array<double, 3> array;
  typedef array::index index;

  array weights;

  zmatrix(const std::vector<Eigen::Vector3d> &, const int);

};

std::vector<double> lagrange_coefficients(int, int, double);
