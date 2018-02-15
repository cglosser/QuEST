#ifndef AIM_INTERACTION_H
#define AIM_INTERACTION_H

#include "interactions/AIM/direct.h"
#include "interactions/AIM/farfield.h"
#include "interactions/AIM/nearfield.h"

namespace AIM {
  class Interaction;
}

class AIM::Interaction final : public InteractionBase {
 public:
  Interaction(const std::shared_ptr<DotVector> dots,
              const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>
                  history,
              Propagation::Kernel<cmplx> &kernel,
              const Eigen::Array3d &spacing,
              const int interp_order,
              const int expansion_order,
              const int border,
              const double c0,
              const double dt,
              normalization::SpatialNorm normalization)
      : InteractionBase(dots, dt),
        grid(spacing, expansion_order, *dots),
        expansion_table(
            AIM::Expansions::LeastSquaresExpansionSolver(expansion_order, grid)
                .table(*dots)),

        ff(dots,
           history,
           interp_order,
           c0,
           dt,
           grid,
           expansion_table,
           normalization),

        nf(dots,
           history,
           interp_order,
           border,
           c0,
           dt,
           grid,
           expansion_table,
           normalization),

        direct(dots, history, kernel, interp_order, border, c0, dt, grid)
  {
  }

  const ResultArray &evaluate(const int t)
  {
    // I DON'T KNOW WHY THAT NEEDS A CONJUGATE!!!
    results =
        (ff.evaluate(t).conjugate() - nf.evaluate(t)) + direct.evaluate(t);
    return results;
  }

 private:
  Grid grid;
  Expansions::ExpansionTable expansion_table;

  Farfield ff;
  Nearfield nf;
  DirectInteraction direct;
};

#endif
