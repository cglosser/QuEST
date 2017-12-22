#ifndef AIM_INTERACTION_H
#define AIM_INTERACTION_H

#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>

#include "common.h"
#include "expansion.h"
#include "fourier.h"
#include "grid.h"
#include "interactions/history_interaction.h"
#include "spacetime.h"

namespace AIM {
  class AimInteraction;

  namespace normalization {
    using SpatialNorm = std::function<double(const Eigen::Vector3d &)>;
    const SpatialNorm unit = [](__attribute__((unused))
                                const Eigen::Vector3d &v) { return 1; };
    const SpatialNorm distance = [](const Eigen::Vector3d &v) {
      return v.norm();
    };
    const SpatialNorm poisson = [](const Eigen::Vector3d &v) {
      return 4 * M_PI * v.norm();
    };
  }
}

class AIM::AimInteraction final : public HistoryInteraction {
 public:
  AimInteraction(const int, const Grid &, normalization::SpatialNorm);
  AimInteraction(
      const std::shared_ptr<const DotVector>,
      const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>,
      const std::shared_ptr<Propagation::RotatingFramePropagator>,
      const int,
      const double,
      const double,
      const Grid &,
      const Expansions::ExpansionTable &,
      Expansions::ExpansionFunction,
      normalization::SpatialNorm);

  const ResultArray &evaluate(const int) final;

  // private:
  Grid grid;
  Expansions::ExpansionTable expansion_table;
  Expansions::ExpansionFunction expansion_function;
  normalization::SpatialNorm normalization;
  std::array<int, 4> circulant_dimensions;

  // This corresponds to delta(t - R/c)/R and thus holds *scalar* quantities
  spacetime::vector<cmplx> fourier_table;

  // These correspond to J and E and thus hold *vector* quantities
  spacetime::vector3d<cmplx> source_table, obs_table;

  TransformPair spatial_vector_transforms;

  void fill_source_table(const int);
  void propagate(const int);
  void fill_results_table(const int);

  spacetime::vector<cmplx> circulant_fourier_table();
  void fill_gmatrix_table(spacetime::vector<cmplx> &) const;
  TransformPair spatial_fft_plans();
};

#endif
