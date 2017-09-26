#include <Eigen/Dense>
#include "../../quantum_dot.h"

namespace AIM {
  class Grid;
}

class AIM::Grid {
 public:
  using BoundsArray = Eigen::Array<int, 3, 2>;
  using ipair_t = std::pair<int, int>;

  Grid();
  Grid(const Eigen::Array3d &, const std::shared_ptr<DotVector> &, const int);
  Grid(const Eigen::Array3d &, const Eigen::Array3i &);

  std::vector<DotRange> box_contents_map(
      const std::shared_ptr<DotVector> &) const;

  // Geometry routines (grid <---> space)
  Eigen::Vector3i grid_coordinate(const Eigen::Vector3d &) const;
  size_t coord_to_idx(const Eigen::Vector3i &) const;
  Eigen::Vector3i idx_to_coord(size_t) const;
  Eigen::Vector3d spatial_coord_of_box(const size_t) const;
  std::vector<size_t> expansion_box_indices(const Eigen::Vector3d &,
                                            const int) const;

  Eigen::Array3i dimensions;
  size_t num_boxes;
  double max_diagonal;

  int max_transit_steps(double c, double dt) const
  {
    return static_cast<int>(ceil(max_diagonal / (c * dt)));
  };

  BoundsArray calculate_bounds() const;

 private:
  Eigen::Array3d spacing;
  std::shared_ptr<DotVector> dots;
  int expansion_order;
  BoundsArray bounds;

  void sort_points_on_boxidx() const;
};
