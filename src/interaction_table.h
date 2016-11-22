#ifndef INTERACTION_TABLE_H
#define INTERACTION_TABLE_H

#include <Eigen/Dense>
#include <algorithm>
#include <boost/multi_array.hpp>
#include <cmath>
#include <vector>

#include "configuration.h"
#include "lagrange_set.h"
#include "quantum_dot.h"

class InteractionTable {
 public:
  InteractionTable(const int, const std::vector<QuantumDot> &);

  int interp_order;
  int num_interactions;

  std::vector<int> floor_delays;
  boost::multi_array<double, 2> coefficients;

  static int coord2idx(int, int);
  static std::pair<int, int> idx2coord(const int);
};

#endif
