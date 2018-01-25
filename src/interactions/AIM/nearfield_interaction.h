#pragma once

#include "interactions/history_interaction.h"
#include "interactions/AIM/grid.h"

namespace AIM {
  class NearfieldInteraction;
}

class AIM::NearfieldInteraction final : public HistoryInteraction {
 public:
  NearfieldInteraction(
      std::shared_ptr<const DotVector>,
      std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>,
      Propagation::RotatingFramePropagator,
      const int,
      const double,
      const double,
      Grid);

  const ResultArray &evaluate(const int) final {
    results.setZero();
    return results;
  };

 //private:
  struct InteractionPair {
    int src_idx, obs_idx;
    std::pair<int, double> delay;
  };

  Propagation::RotatingFramePropagator propagator;
  int num_interactions;
  std::vector<InteractionPair> interaction_pairs;
  boost::multi_array<cmplx, 2> coefficients;
  Grid grid;

  void build_pair_list();
};

