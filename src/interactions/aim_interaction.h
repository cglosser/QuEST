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
  typedef std::pair<DotVector::iterator, DotVector::iterator> BoxRange;

  Grid(const Eigen::Vector3d &, const std::shared_ptr<DotVector> &);

  Eigen::Vector3i num_boxes;
  std::vector<BoxRange> boxes;

 private:
  Eigen::Vector3d spacing;
  std::shared_ptr<DotVector> dots;
  BoundsArray bounds;

  BoundsArray calculate_bounds() const;
  void sort_points_on_boxidx() const;
  void map_points_to_boxes();
  Eigen::Vector3i grid_coordinate(const Eigen::Vector3d &) const;
  size_t coord_to_idx(const Eigen::Vector3i &) const;
  Eigen::Vector3i idx_to_coord(size_t) const;
};
#endif
