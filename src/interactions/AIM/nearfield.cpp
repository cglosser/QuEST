#include "nearfield.h"

AIM::Nearfield::Nearfield(
    const std::shared_ptr<const DotVector> dots,
    const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history,
    const int interp_order,
    const int border,
    const double c0,
    const double dt,
    const Grid grid,
    const Expansions::ExpansionTable &expansion_table,
    Expansions::ExpansionFunction expansion_function,
    normalization::SpatialNorm normalization)
    : AimBase(dots,
              history,
              interp_order,
              c0,
              dt,
              grid,
              expansion_table,
              expansion_function,
              normalization,
              grid.nearfield_shape(c0, dt, interp_order, border)),
    neighbors(grid.nearfield_pairs(border))
{
  propagation_table = make_propagation_table();
}
