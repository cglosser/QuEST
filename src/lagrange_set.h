#ifndef LAGRANGE_SET_H
#define LAGRANGE_SET_H

#include <boost/multi_array.hpp>

class UniformLagrangeSet
{
 public:
  typedef boost::multi_array<double, 2> Array;

  int order;
  Array weights;

  explicit UniformLagrangeSet(const int);
  UniformLagrangeSet(const int, const double);

  void calculate_weights(const double);
};

#endif
