#ifndef AIM_INTERACTION_H
#define AIM_INTERACTION_H

#include <Eigen/Dense>
#include <algorithm>

#include "../quantum_dot.h"

namespace AIM {
  class Grid;
}

class AIM::Grid {
 public:
  typedef Eigen::Array<int, 3, 2> BoundsArray;

  Grid(const Eigen::Vector3d &, const std::shared_ptr<DotVector> &);

  std::vector<std::vector<size_t>> point_mapping;

 private:
  Eigen::Vector3d spacing;
  std::shared_ptr<DotVector> dots;
  BoundsArray bounds;
  Eigen::Vector3i num_boxes;

  BoundsArray calculate_bounds() const;
  void map_points_to_boxes();
  void sort_points_on_boxidx() const;
  Eigen::Vector3i grid_coordinate(const Eigen::Vector3d &) const;
  size_t coord_to_idx(const Eigen::Vector3i &) const;
  Eigen::Vector3i idx_to_coord(size_t) const;
};
#endif
