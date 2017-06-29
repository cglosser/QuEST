#ifndef HISTORY_INTERACTION_H
#define HISTORY_INTERACTION_H

#include <Eigen/Dense>
#include <boost/multi_array.hpp>

#include "../integrator/history.h"
#include "../lagrange_set.h"
#include "../quantum_dot.h"
#include "green_function.h"
#include "interaction.h"

class HistoryInteraction : public Interaction {
 public:
  HistoryInteraction(
      const std::shared_ptr<const DotVector> &,
      const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> &,
      const std::shared_ptr<Propagation::RotatingFramePropagator> &, const int,
      const double, const double);

  virtual const ResultArray &evaluate(const int);

 private:
  std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history;
  std::shared_ptr<Propagation::RotatingFramePropagator> dyadic;
  int interp_order, num_interactions;
  std::vector<int> floor_delays;
  boost::multi_array<cmplx, 2> coefficients;
  const double dt;
  const double c0;

  void build_coefficient_table();

  static int coord2idx(int, int);
  static std::pair<int, int> idx2coord(const int);
  static std::pair<int, double> split_double(const double);
};

#endif
