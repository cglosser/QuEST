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
      grid(std::move(grid))
{
  build_pair_list();
}

void AIM::NearfieldInteraction::build_pair_list()
{
  auto box_contents = grid.box_contents_map();

  for(auto src_idx = 0u; src_idx < grid.num_gridpoints - 1; ++src_idx) {
    Eigen::Vector3i src_coord = grid.idx_to_coord(src_idx);
    for(auto obs_idx = src_idx + 1; obs_idx < grid.num_gridpoints; ++obs_idx) {
      Eigen::Vector3i obs_coord = grid.idx_to_coord(obs_idx);
      Eigen::Vector3i dr = obs_coord - src_coord;
      bool nearfield_point = dr.minCoeff() < interp_order;
      if(!nearfield_point)
        continue;  // complement of condition in AimInteraction

      for(auto src_dot = box_contents[src_idx].first;
          src_dot != box_contents[src_idx].second; ++src_dot) {

        int p1 = std::distance<DotVector::const_iterator>(dots->begin(), src_dot);

        for(auto obs_dot = box_contents[obs_idx].first;
            obs_dot != box_contents[obs_idx].second; ++obs_dot) {

          int p2 = std::distance<DotVector::const_iterator>(dots->begin(), obs_dot);

          Eigen::Vector3d dr(separation((*dots)[p1], (*dots)[p2]));
          std::pair<int, double> delay(split_double(dr.norm() / (c0 * dt)));

          interaction_pairs.push_back({p1, p2, delay});

        }
      }
    }
  }
}
