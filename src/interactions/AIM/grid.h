#ifndef GRID_H
#define GRID_H

#include <Eigen/Dense>
#include "math_utils.h"
#include "quantum_dot.h"

namespace AIM {
  class Grid;
}

class AIM::Grid {
 public:
  using BoundsArray = Eigen::Array<int, 3, 2>;
  using ipair_t = std::pair<int, int>;

  Grid();
  Grid(const Eigen::Array3d &, const std::shared_ptr<DotVector>, const int);
  Grid(const Eigen::Array3d &,
       const Eigen::Array3i &,
       const Eigen::Vector3i & = Eigen::Vector3i::Zero());

  BoundsArray calculate_bounds() const;
  std::array<int, 4> circulant_shape(const double,
                                     const double,
                                     const int = 0) const;
  std::vector<DotRange> box_contents_map() const;
  int expansion_order() const { return expansion_order_; }
  // Geometry routines (grid <---> space)
  inline Eigen::Vector3i grid_coordinate(const Eigen::Vector3d &coord) const
  {
    return (coord.array() / spacing).cast<int>();
  }

  inline size_t associated_grid_index(const Eigen::Vector3d &coord) const
  {
    Eigen::Vector3i grid_coord = grid_coordinate(coord);
    return coord_to_idx(grid_coord - bounds.col(0).matrix());
  }

  inline size_t coord_to_idx(const Eigen::Vector3i &coord) const
  {
    return coord(2) + dimensions(2) * (coord(1) + dimensions(1) * coord(0));
  }

  inline Eigen::Vector3i idx_to_coord(size_t idx) const
  {
    const int nynz = dimensions(1) * dimensions(2);
    const int x = idx / nynz;
    idx -= x * nynz;
    const int y = idx / dimensions(2);
    const int z = idx % dimensions(2);

    return Eigen::Vector3i(x, y, z);
  }

  inline Eigen::Vector3d spatial_coord_of_box(const size_t box_id) const
  {
    Eigen::Vector3i dr = (idx_to_coord(box_id) + bounds.col(0).matrix());
    return dr.array().cast<double>() * spacing;
  }
  std::vector<size_t> expansion_box_indices(const Eigen::Vector3d &) const;

  int min_distance(const Eigen::Vector3d &r1, const Eigen::Vector3d &r2) const
  {
    auto r1_grid_indices = expansion_box_indices(r1);
    auto r2_grid_indices = expansion_box_indices(r2);

    int min_dist =
        (idx_to_coord(r1_grid_indices[0]) - idx_to_coord(r2_grid_indices[0]))
            .lpNorm<Eigen::Infinity>();
    for(const auto &r1_idx : r1_grid_indices) {
      for(const auto &r2_idx : r2_grid_indices) {
        Eigen::Vector3i dr = idx_to_coord(r1_idx) - idx_to_coord(r2_idx);
        min_dist = std::min(min_dist, dr.lpNorm<Eigen::Infinity>());
        if(min_dist == 0) return min_dist;  // Not going to find a smaller one
      }
    }

    return min_dist;
  }

  bool is_nearfield_pair(const Eigen::Vector3d &r1,
                         const Eigen::Vector3d &r2) const
  {
    return min_distance(r1, r2) < expansion_order_;
  }

  Eigen::Array3i dimensions;
  size_t num_gridpoints;
  double max_diagonal;

  int max_transit_steps(double c, double dt) const
  {
    return static_cast<int>(ceil(max_diagonal / (c * dt)));
  };

 private:
  Eigen::Array3d spacing;
  std::shared_ptr<DotVector> dots;
  int expansion_order_;
  BoundsArray bounds;

  void sort_points_on_boxidx() const;
};

#endif
