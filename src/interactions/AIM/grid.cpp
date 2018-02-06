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

std::array<int, 4> AIM::Grid::toeplitz_shape(const double c,
                                             const double dt,
                                             const int pad /* = 0 */) const
{
  std::array<int, 4> shape;

  shape[0] = max_transit_steps(c, dt) + pad;
  for(int i = 1; i < 4; ++i) shape[i] = dimensions(i - 1);

  return shape;
}

std::array<int, 4> AIM::Grid::circulant_shape(const double c,
                                              const double dt,
                                              const int pad /* = 0 */) const
{
  std::array<int, 4> shape(toeplitz_shape(c, dt, pad));
  for(int i = 1; i < 4; ++i) shape[i] *= 2;

  return shape;
}

std::vector<DotRange> AIM::Grid::box_contents_map() const
{
  std::vector<DotRange> boxes(num_gridpoints);
  for(size_t box_idx = 0; box_idx < boxes.size(); ++box_idx) {
    auto NearestGridpoint = [=](const QuantumDot &qd) {
      return associated_grid_index(qd.position()) == box_idx;
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

std::set<size_t> AIM::Grid::active_nodes() const
{
  std::set<size_t> active;

  const auto mapping = box_contents_map();
  for(int i = 0; i < static_cast<int>(mapping.size()); ++i) {
    if(mapping[i].first == mapping[i].second) continue;

    const auto expansion_pts = expansion_indices(i);
    active.insert(expansion_pts.begin(), expansion_pts.end());
  }

  return active;
}

std::vector<AIM::Grid::ipair_t> AIM::Grid::full_nearfield_pairs(
    const int thresh) const
{
  // List of gridpoints...
  //    ...that have an expansion...
  //        within thresh of each other. Thresh is in grid units.

  std::vector<ipair_t> nf;
  if(!dots) return nf;

  const int ubound = expansion_order + thresh;
  const int lbound = -ubound;
  const auto mapping = box_contents_map();
  const auto active = active_nodes();

  for(const auto &src_idx : active) {
    Eigen::Vector3i src_coord = idx_to_coord(src_idx);
    for(int x = lbound; x <= ubound; ++x) {
      for(int y = lbound; y <= ubound; ++y) {
        for(int z = lbound; z <= ubound; ++z) {
          const Eigen::Vector3i delta(x, y, z);
          const Eigen::Vector3i new_pt = src_coord + delta;
          const auto obs_idx = coord_to_idx(new_pt);

          if((new_pt.array() < 0).any() ||
             (new_pt.array() >= dimensions).any() || obs_idx <= src_idx ||
             (!active.count(obs_idx)))
            continue;

          nf.emplace_back(src_idx, obs_idx);
        }
      }
    }
  }

  nf.shrink_to_fit();
  return nf;
}

std::vector<AIM::Grid::ipair_t> AIM::Grid::compressed_nearfield_pairs(
    const int thresh) const
{
  // Same as full_nearfield_pairs except both gridpoints must also be
  // "active" (have associated particles).

  std::vector<ipair_t> nf;
  if(!dots) return nf;

  const int ubound = expansion_order + thresh;
  const int lbound = -ubound;
  const auto mapping = box_contents_map();

  for(auto src_idx = 0u; src_idx < num_gridpoints; ++src_idx) {
    if(mapping[src_idx].first == mapping[src_idx].second) continue;

    Eigen::Vector3i src_coord = idx_to_coord(src_idx);
    for(int x = lbound; x <= ubound; ++x) {
      for(int y = lbound; y <= ubound; ++y) {
        for(int z = lbound; z <= ubound; ++z) {
          const Eigen::Vector3i delta(x, y, z);
          const Eigen::Vector3i new_pt = src_coord + delta;
          const auto obs_idx = coord_to_idx(new_pt);

          const bool invalid_point = (new_pt.array() < 0).any() ||
                                     (new_pt.array() >= dimensions).any();

          const bool empty_obs =
              mapping[obs_idx].first == mapping[obs_idx].second;

          if(invalid_point || empty_obs || obs_idx < src_idx) continue;

          nf.emplace_back(src_idx, obs_idx);
        }
      }
    }
  }

  nf.shrink_to_fit();
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
