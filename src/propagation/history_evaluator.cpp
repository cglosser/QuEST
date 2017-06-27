#include "history_evaluator.h"

template <class scalarFieldType, class soltype>
HistoryEvaluator::HistoryEvaluator(
    const std::shared_ptr<const Integrator::History<soltype>> &history,
    const std::shared_ptr<GreenFunction::Dyadic> &dyad, const int interp_order,
    const double c0, const double dt) :
      FieldEvaluator<scalarFieldType>(dots), history(history) 
{}
