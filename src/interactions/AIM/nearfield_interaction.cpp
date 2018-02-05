#include "nearfield_interaction.h"

AIM::NearfieldInteraction::NearfieldInteraction(
    std::shared_ptr<const DotVector> dots,
    std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history,
    Propagation::RotatingFramePropagator propagator,
    const int interp_order,
    const double c0,
    const double dt,
    Grid grid)
    : HistoryInteraction(
          std::move(dots), std::move(history), interp_order, c0, dt),
      propagator(std::move(propagator)),
      grid(std::move(grid)),
      floor_delays(interaction_pairs.size()),
      coefficients(boost::extents[interaction_pairs.size()][interp_order + 1])
{
  build_coefficient_table();
}

void AIM::NearfieldInteraction::build_coefficient_table()
{
  Interpolation::UniformLagrangeSet lagrange(interp_order);

  for(auto pair_idx = 0u; pair_idx < interaction_pairs.size(); ++pair_idx) {
    std::pair<int, int> &pair = interaction_pairs[pair_idx];

    Eigen::Vector3d dr(separation((*dots)[pair.first], (*dots)[pair.second]));
    auto delay = split_double(dr.norm() / (c0 * dt));

    floor_delays[pair_idx] = delay.first;
    lagrange.evaluate_derivative_table_at_x(delay.second, dt);

    std::vector<Eigen::Matrix3cd> interp_dyads(
        propagator.coefficients(dr, lagrange));

    for(int i = 0; i <= interp_order; ++i) {
      coefficients[pair_idx][i] = (*dots)[pair.second].dipole().dot(
          interp_dyads[i] * (*dots)[pair.first].dipole());
    }
  }
}

const Interaction::ResultArray &AIM::NearfieldInteraction::evaluate(
    const int time_idx)
{
  results.setZero();

  for(auto pair_idx = 0u; pair_idx < interaction_pairs.size(); ++pair_idx) {
    const auto &pair = interaction_pairs[pair_idx];
    const int s = time_idx - floor_delays[pair_idx];

    for(int i = 0; i <= interp_order; ++i) {
      if(s - i < history->array.index_bases()[1]) continue;

      constexpr int RHO_01 = 1;

      results[pair.first] += (history->array[pair.second][s - i][0])[RHO_01] *
                             coefficients[pair_idx][i];
      results[pair.second] += (history->array[pair.first][s - i][0])[RHO_01] *
                              coefficients[pair_idx][i];
    }
  }

  return results;
}
