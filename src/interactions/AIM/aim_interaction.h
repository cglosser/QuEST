#ifndef AIM_INTERACTION_H
#define AIM_INTERACTION_H

#include <fftw3.h>
#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>

#include <iomanip>

#include "common.h"
#include "expansion.h"
#include "grid.h"
#include "interactions/history_interaction.h"

namespace AIM {
  class AimInteraction;
}

class AIM::AimInteraction final : public HistoryInteraction {
 public:

  struct TransformPair {
    fftw_plan forward, backward;
    ~TransformPair()
    {
      fftw_destroy_plan(forward);
      fftw_destroy_plan(backward);
    }
  };

  AimInteraction(
      const std::shared_ptr<const DotVector> &,
      const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> &,
      const std::shared_ptr<Propagation::RotatingFramePropagator> &,
      const int,
      const double,
      const double,
      const Grid &,
      const Array<Expansion> &);

  const ResultArray &evaluate(const int) final;


  // private:
  Grid grid;
  Array<Expansion> expansion_table;
  int box_order, max_transit_steps;
  std::array<int, 4> circulant_dimensions;

  SpacetimeVector<cmplx> fourier_table, source_table, obs_table;
  TransformPair spatial_transforms;

  void fill_source_table(const int);

  SpacetimeVector<cmplx> circulant_fourier_table();
  void fill_gmatrix_table(SpacetimeVector<cmplx> &) const;
  TransformPair spatial_fft_plans();
};

#endif
