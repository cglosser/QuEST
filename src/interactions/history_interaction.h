#ifndef HISTORY_INTERACTION_H
#define HISTORY_INTERACTION_H

#include <Eigen/Dense>
#include <boost/multi_array.hpp>

#include "../integrator/history.h"
#include "../lagrange_set.h"
#include "../math_utils.h"
#include "../quantum_dot.h"
#include "green_function.h"
#include "interaction.h"

class HistoryInteraction : public Interaction {
 public:
  HistoryInteraction(
      const std::shared_ptr<const DotVector> &dots,
      const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>
          &history,
      const std::shared_ptr<Propagation::RotatingFramePropagator> &propagator,
      const int interp_order,
      const double c0,
      const double dt)
      : Interaction(dots, dt),
        history(history),
        propagator(propagator),
        interp_order(interp_order),
        c0(c0){};

 protected:
  std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history;
  std::shared_ptr<Propagation::RotatingFramePropagator> propagator;
  int interp_order;
  const double c0;
};

#endif
