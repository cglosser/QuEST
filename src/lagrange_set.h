#ifndef LAGRANGE_SET_H
#define LAGRANGE_SET_H

#include <boost/multi_array.hpp>

namespace Interpolation {
  class UniformLagrangeSet;
  class DerivFive;
  typedef boost::multi_array<double, 2> InterpolationTable;
  constexpr int NUM_DERIVATIVES = 3;
}

class Interpolation::UniformLagrangeSet {
 public:
  InterpolationTable evaluations;

  UniformLagrangeSet(const int);
  UniformLagrangeSet(const double, const int, const double = 1);

  void evaluate_derivative_table_at_x(const double, const double = 1);
  int order() const { return order_; }
 private:
  int order_;
};

class Interpolation::DerivFive {
 public:
  InterpolationTable evaluations;

  DerivFive(const double = 1);

  void evaluate_derivative_table_at_x(const double);

  double d0(const int, const double) const;
  double d1(const int, const double) const;
  double d2(const int, const double) const;

  constexpr static int order() { return 5; }
 private:
  double dt_;
};

#endif
