#include "grid.h"

AIM::Grid::Grid()
    : Grid(Eigen::Vector3d(0, 0, 0), std::make_shared<DotVector>(), 0)
{
}

AIM::Grid::Grid(const Eigen::Array3d &spacing,
                const std::shared_ptr<DotVector> &dots,
                const int expansion_order)
    : spacing(spacing),
      dots(dots),
      expansion_order(expansion_order),
      bounds(calculate_bounds())
{
  dimensions = bounds.col(1) - bounds.col(0) + 1;
  num_boxes = dimensions.prod();
  max_diagonal = (dimensions.cast<double>() * spacing).matrix().norm();

  sort_points_on_boxidx();
}

AIM::Grid::Grid(const Eigen::Array3d &spacing, const Eigen::Array3i &dimensions)
    : dimensions(dimensions),
      spacing(spacing),
      dots(nullptr),
      expansion_order(0)
{
  bounds.col(0) = 0;
  bounds.col(1) = dimensions;

  num_boxes = dimensions.prod();
  max_diagonal = (dimensions.cast<double>() * spacing).matrix().norm();
}

AIM::Grid::BoundsArray AIM::Grid::calculate_bounds() const
{
  BoundsArray b;
  b.col(0) = std::numeric_limits<int>::max();
  b.col(1) = std::numeric_limits<int>::min();

  for(const auto &qdot : *dots) {
    Eigen::Vector3i grid_coord = grid_coordinate(qdot.position());

    b.col(0) = grid_coord.array().min(b.col(0));
    b.col(1) = grid_coord.array().max(b.col(1));
  }

  b.col(0) -= expansion_order / 2;
  b.col(1) += (expansion_order + 1) / 2;

  return b;
}

std::pair<AIM::Grid::ipair_t, AIM::Grid::ipair_t> AIM::Grid::fourier_idx_pairs(
    const int row, const int col) const
{
  auto src = std::make_pair(col / dimensions(1), col % dimensions(1));
  auto obs = std::make_pair(row / dimensions(1), row % dimensions(1));

  return std::make_pair(src, obs);
}

int AIM::Grid::fourier_idx(const int row, const int col) const
{
  ipair_t src, obs;
  std::tie(src, obs) = fourier_idx_pairs(row, col);
  return dimensions(1) * std::abs(obs.first - src.first) +
         std::abs(obs.second - src.second);
}
