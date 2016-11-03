#ifndef LAGRANGE_SET_H
#define LAGRANGE_SET_H

#include <boost/multi_array.hpp>

#include "configuration.h"

class UniformLagrangeSet
{
  public:
  UniformLagrangeSet(const double);

  double sample_x;
  boost::multi_array<double, 2> weights;

  private:
  void calculate_weights();
};

#endif
