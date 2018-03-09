#ifndef AIM_BASE_H
#define AIM_BASE_H

#include "chebyshev.h"
#include "expansion.h"
#include "grid.h"
#include "interactions/history_interaction.h"
#include "normalization.h"
#include "projector.h"
#include "spacetime.h"

namespace AIM {
  class AimBase;

  constexpr int chebyshev_order = 4;
}

class AIM::AimBase : public HistoryInteraction {
 public:
  AimBase(const std::shared_ptr<const DotVector> dots,
          const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>
              history,
          const int interp_order,
          const double c0,
          const double dt,
          const Grid grid,
          const Expansions::ExpansionTable expansion_table,
          Normalization::SpatialNorm normalization,
          std::array<int, 4> table_dimensions,
          const boost::multi_array<double, 4> chebyshev_weights,
          Projector::Projector_fn<cmplx> projector)
      : HistoryInteraction(dots, history, interp_order, c0, dt),

        grid(std::move(grid)),
        expansion_table(std::move(expansion_table)),
        normalization(std::move(normalization)),
        table_dimensions(table_dimensions),

        // Tables
        propagation_table(table_dimensions),
        source_table(spacetime::make_vector3d<cmplx>(table_dimensions)),
        obs_table(spacetime::make_vector3d<cmplx>(table_dimensions)),

        // Interior grid
        Cheb(grid, *dots),
        projector_(std::move(projector)),
        chebyshev_weights(chebyshev_weights),
        chebyshev_table(boost::extents[table_dimensions[0]][grid.size()]
                                      [chebyshev_order + 1][chebyshev_order + 1]
                                      [chebyshev_order + 1][3]),
        chebyshev_coefficients(
            boost::extents[table_dimensions[0]][grid.size()]
                          [chebyshev_order + 1][chebyshev_order + 1]
                          [chebyshev_order + 1][3])
  {
    auto clear = [](auto &table) {
      std::fill(table.data(), table.data() + table.num_elements(), cmplx(0, 0));
    };

    clear(propagation_table);
    clear(source_table);
    clear(obs_table);
  }

  const ResultArray &evaluate(const int step) final
  {
    fill_source_table(step);
    propagate(step);
    fill_chebyshev_table(step);

    const auto &x = Cheb.interpolate(step, chebyshev_coefficients, projector_);

    for(size_t i = 0; i < dots->size(); ++i) {
      results(i) = x.col(i).dot((*dots)[i].dipole());
    }

    return results;
  }

 protected:
  const Grid grid;
  const Expansions::ExpansionTable expansion_table;
  Normalization::SpatialNorm normalization;
  boost::multi_array<int, 2> expansion_indices;

  std::array<int, 4> table_dimensions;

  // This corresponds to delta(t - R/c)/R and thus holds *scalar*
  // quantities
  spacetime::vector<cmplx> propagation_table;

  // These correspond to J and E and thus hold *vector* quantities
  spacetime::vector3d<cmplx> source_table, obs_table;

  // Interior grid
  Chebyshev<cmplx, chebyshev_order> Cheb;
  Projector::Projector_fn<cmplx> projector_;
  const boost::multi_array<double, 4> chebyshev_weights;
  boost::multi_array<cmplx, 6> chebyshev_table, chebyshev_coefficients;

  virtual spacetime::vector<cmplx> make_propagation_table() const = 0;
  virtual void fill_source_table(const int) = 0;
  virtual void propagate(const int) = 0;
  virtual void fill_chebyshev_table(const int) = 0;
};

#endif
