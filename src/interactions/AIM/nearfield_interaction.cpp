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
      interaction_pairs(build_pair_list()),
      floor_delays(interaction_pairs.size()),
      coefficients(boost::extents[interaction_pairs.size()][interp_order + 1])
{
  build_coefficient_table();
}

std::vector<std::pair<int, int>> AIM::NearfieldInteraction::build_pair_list()
    const
{
  std::vector<std::pair<int, int>> pairs;
  auto box_contents = grid.box_contents_map();

  for(auto i = 0u; i < box_contents.size() - 1; ++i) {
    // start iter = end iter, thus empty box
    if(box_contents[i].first == box_contents[i].second) continue;

    for(auto j = i; j < box_contents.size(); ++j) {
      if(box_contents[j].first == box_contents[j].second)
        continue;  // start iter = end iter, thus empty box

      // box_contents[i] and [j] yield *pairs of DotVector iterators*
      // corresponding to the range of particles within the box (which assumes
      // the DotVector is sorted). If the *.first and *.second iterators are
      // equal, then the box is empty (checked above), otherwise it contains
      // particles that can ALL equivalently determine the box's position and
      // thus its nearfield neighbors.
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
            pairs.push_back({p1, p2});
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
            pairs.push_back({p1, p2});
          }
        }
      }
    }
  }

  pairs.shrink_to_fit();
  return pairs;
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
