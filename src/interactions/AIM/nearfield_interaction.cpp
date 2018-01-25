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

  for(auto i = 0u; i < box_contents.size() - 1; ++i) {
    // start iter = end iter, thus empty box
    if(box_contents[i].first == box_contents[i].second) continue;

    for(auto j = i; j < box_contents.size(); ++j) {
      if(box_contents[j].first == box_contents[j].second)
        continue;  // start iter = end iter, thus empty box

      // box_contents[i,j] have particles, so the first particle in
      // each can be used to determine box separation/nearfieldness
      bool is_in_nearfield = grid.is_nearfield_pair(
          box_contents[i].first->position(), box_contents[j].first->position());
      if(!is_in_nearfield) continue;

      if(i == j) {
        // i and j point to the same box; this looping
        // structure avoids double counting pairs

        for(auto src_dot = box_contents[i].first;
            src_dot != box_contents[i].second - 1; ++src_dot) {
          for(auto obs_dot = src_dot + 1; obs_dot != box_contents[i].second;
              ++obs_dot) {
            int p1 = std::distance<DotVector::const_iterator>(dots->begin(),
                                                              src_dot);
            int p2 = std::distance<DotVector::const_iterator>(dots->begin(),
                                                              obs_dot);

            Eigen::Vector3d dr(separation((*dots)[p1], (*dots)[p2]));
            auto delay = split_double(dr.norm() / (c0 * dt));

            interaction_pairs.push_back({p1, p2, delay});
          }
        }

      } else {
        // i and j are different boxes, ergo we need
        // to loop over all particles in both
        for(auto src_dot = box_contents[i].first;
            src_dot != box_contents[i].second; ++src_dot) {
          for(auto obs_dot = box_contents[j].first;
              obs_dot != box_contents[j].second; ++obs_dot) {
            int p1 = std::distance<DotVector::const_iterator>(dots->begin(),
                                                              src_dot);
            int p2 = std::distance<DotVector::const_iterator>(dots->begin(),
                                                              obs_dot);

            Eigen::Vector3d dr(separation((*dots)[p1], (*dots)[p2]));
            auto delay = split_double(dr.norm() / (c0 * dt));

            interaction_pairs.push_back({p1, p2, delay});
          }
        }
      }
    }
  }
}
