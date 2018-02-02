#include "grid.h"

AIM::Grid::Grid() : Grid(Eigen::Vector3d(0, 0, 0), nullptr, 0) {}
AIM::Grid::Grid(const Eigen::Array3d &spacing,
                const std::shared_ptr<DotVector> dots,
                const int expansion_order)
    : spacing(spacing),
      dots(std::move(dots)),
      expansion_order(expansion_order),
      bounds(calculate_bounds())
{
  dimensions = bounds.col(1) - bounds.col(0) + 1;
  num_gridpoints = dimensions.prod();
  max_diagonal = ((bounds.col(1) - bounds.col(0)).cast<double>() * spacing)
                     .matrix()
                     .norm();

  sort_points_on_boxidx();
}

AIM::Grid::Grid(const Eigen::Array3d &spacing,
                const Eigen::Array3i &dimensions,
                const Eigen::Vector3i &shift,
                const int expansion_order)
    : dimensions(dimensions),
      spacing(spacing),
      dots(nullptr),
      expansion_order(expansion_order)
{
  bounds.col(0) = 0 + shift.array();
  bounds.col(1) = dimensions + shift.array() - 1;

  num_gridpoints = dimensions.prod();
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

std::array<int, 4> AIM::Grid::circulant_shape(const double c,
                                              const double dt,
                                              const int pad /* = 0 */) const
{
  std::array<int, 4> shape;

  shape[0] = max_transit_steps(c, dt) + pad;
  for(int i = 1; i < 4; ++i) shape[i] = 2 * dimensions(i - 1);

  return shape;
}

std::vector<DotRange> AIM::Grid::box_contents_map(
    const std::shared_ptr<DotVector> &dots) const
{
  std::vector<DotRange> boxes(num_gridpoints);
  for(size_t box_idx = 0; box_idx < boxes.size(); ++box_idx) {
    auto NearestGridpoint = [=](const QuantumDot &qd) {
      return coord_to_idx(grid_coordinate(qd.position())) == box_idx;
    };

    auto begin = std::find_if(dots->begin(), dots->end(), NearestGridpoint);
    auto end = std::find_if_not(begin, dots->end(), NearestGridpoint);
    boxes.at(box_idx) = std::make_pair(begin, end);
  }

  return boxes;
}

std::vector<size_t> AIM::Grid::expansion_indices(const int grid_index) const
{
  Eigen::Vector3i origin = idx_to_coord(grid_index);
  std::vector<size_t> indices(std::pow(expansion_order + 1, 3));

  size_t idx = 0;
  for(int nx = 0; nx <= expansion_order; ++nx) {
    for(int ny = 0; ny <= expansion_order; ++ny) {
      for(int nz = 0; nz <= expansion_order; ++nz) {
        const Eigen::Vector3i delta(grid_sequence(nx), grid_sequence(ny),
                                    grid_sequence(nz));
        const size_t grid_idx = coord_to_idx(origin + delta);

        indices.at(idx++) = grid_idx;
      }
    }
  }

  return indices;
}

int AIM::Grid::expansion_distance(const int i, const int j) const
{
  auto pts1 = expansion_indices(i);
  auto pts2 = expansion_indices(j);

  int min_dist = (idx_to_coord(pts1.front()) - idx_to_coord(pts2.front()))
                     .lpNorm<Eigen::Infinity>();
  for(const auto &p1 : pts1) {
    for(const auto &p2 : pts2) {
      Eigen::Vector3i dr = idx_to_coord(p2) - idx_to_coord(p1);
      min_dist = std::min(min_dist, dr.lpNorm<Eigen::Infinity>());
      if(min_dist == 0) return 0;
    }
  }

  return min_dist;
}

std::vector<AIM::Grid::ipair_t> AIM::Grid::nearfield_pairs(
    const int thresh) const
{
  // Threshold is in grid units

  std::vector<ipair_t> nf;
  if(!dots) return nf;

  const auto mapping = box_contents_map(dots);
  for(auto pt1 = 0u; pt1 < num_gridpoints - 1; ++pt1) {
    if(mapping[pt1].first == mapping[pt1].second) continue;  // empty box
    auto r1 = idx_to_coord(pt1);

    for(auto pt2 = pt1; pt2 < num_gridpoints; ++pt2) {
      if(mapping[pt2].first == mapping[pt2].second) continue;
      auto r2 = idx_to_coord(pt2);

      auto dist = (r2 - r1).lpNorm<Eigen::Infinity>();
      if(dist < thresh) {
        nf.push_back({pt1, pt2});
      }
    }
  }

  return nf;
}

void AIM::Grid::sort_points_on_boxidx() const
{
  auto grid_comparitor = [&](const QuantumDot &q1, const QuantumDot &q2) {
    return associated_grid_index(q1.position()) <
           associated_grid_index(q2.position());
  };

  std::stable_sort(dots->begin(), dots->end(), grid_comparitor);
}
