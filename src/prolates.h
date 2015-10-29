#ifndef PROLATES_H
#define PROLATES_H

#include <Eigen/Dense>
#include <Eigen/StdVector> // Special allocator for std::vectors
#include <boost/math/constants/constants.hpp>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <iostream>

constexpr double TOLER = 1e-5;

class Prolate {
 private:
  double alpha;
  int width;

  double sqrt_term(double) const;
  double d1_sqrt_term(double) const;
  double d2_sqrt_term(double) const;
  double d3_sqrt_term(double) const;

 public:
  Prolate(int);

  int get_width() {return width;}

  double d0(double) const; //0th derivative
  double d1(double) const; //1st derivative
  double d2(double) const; //2nd derivative
};

class ProlateTimeExpansion {
 public:
  ProlateTimeExpansion(const Prolate &, const int, const Eigen::Vector3d &);
  void step(const Eigen::Vector3d &);
  Eigen::Vector3d at(const double);

 private:
  Prolate basis_function;
  std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d>> history;


};

#endif
