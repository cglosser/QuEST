#include <boost/multi_array.hpp>

#include "configuration.h"

class UniformLagrangeSet
{
  public:
  UniformLagrangeSet(const double);

  double sample_x;
  boost::multi_array<double, 2> weights;

  private:
  void set_deriv_0();
  void set_deriv_1();
  void set_deriv_2();
};
