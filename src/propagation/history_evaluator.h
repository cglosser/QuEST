#ifndef HISTORY_EVALUATOR_H
#define HISTORY_EVALUATOR_H

#include "propagation.h"
#include "../integrator/history.h"
#include "green_function.h"

namespace Propagation {
  template <class scalarFieldType, class soltype>
  class HistoryEvaluator;
}

template <class scalarFieldType, class soltype>
class Propagation::HistoryEvaluator
    : public Propagation::FieldEvaluator<scalarFieldType> {
 public:
  HistoryEvaluator(const std::shared_ptr<const DotVector> &,
                   const std::shared_ptr<const Integrator::History<soltype>> &,
                   const std::shared_ptr<GreenFunction::Dyadic> &, const int,
                   const double, const double);

 private:
  std::shared_ptr<const Integrator::History<soltype>> history;
  std::shared_ptr<GreenFunction::Dyadic> dyadic;
  int interp_order, num_interactions;
  std::vector<int> floor_delays;
  boost::multi_array<Eigen::Matrix<scalarFieldType, 3, 3>, 2>
      interpolated_propagators;
  const double c0, dt;

  void build_coefficient_table();

  static int coord2idx(const std::pair<int, int> &);
  static std::pair<int, int> idx2coord(const int);
  static std::pair<int, double> split_double(const double);
};

#endif
