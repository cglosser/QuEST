#include <iostream>
#include <Eigen/Dense>
#include <boost/multi_array.hpp>

class UniformLagrangeSet
{
  public:
  UniformLagrangeSet(const double, const int);

  double sample_x;
  int order;
  boost::multi_array<double, 2> coefficients;

  private:
  void set_deriv_0();
  void set_deriv_1();
  void set_deriv_2();
};
