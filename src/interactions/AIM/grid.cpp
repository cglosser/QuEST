#include "grid.h"

AIM::Grid::Grid(const Eigen::Array3d &spacing,
                const int order,
                const Eigen::Array3i &dimensions,
                const Eigen::Vector3i &shift)
    : spacing_(spacing), order_(order), dimensions(dimensions)
{
  bounds.col(0) = 0 + shift.array();
  bounds.col(1) = dimensions + shift.array() - 1;

  num_gridpoints = dimensions.prod();
}

AIM::Grid::Grid(const Eigen::Array3d &spacing, const int order, DotVector &dots)
    : spacing_(spacing),
      order_(order),
      bounds(calculate_bounds(dots)),
      dimensions(bounds.col(1) - bounds.col(0) + 1),
      num_gridpoints(dimensions.prod())
{
  sort_points_on_boxidx(dots);
}

AIM::Grid::BoundsArray AIM::Grid::calculate_bounds(const DotVector &dots) const
{
  BoundsArray b;
  b.col(0) = std::numeric_limits<int>::max();
  b.col(1) = std::numeric_limits<int>::min();

  for(const auto &qdot : dots) {
    Eigen::Vector3i grid_coord = grid_coordinate(qdot.position());

    b.col(0) = grid_coord.array().min(b.col(0));
    b.col(1) = grid_coord.array().max(b.col(1));
  }

  b.col(0) -= order_ / 2;
  b.col(1) += (order_ + 1) / 2;

  return b;
}

std::array<int, 4> AIM::Grid::circulant_shape(const double c,
                                              const double dt,
                                              const int pad /* = 0 */) const
{
  std::array<int, 4> dims;
  dims[0] = max_transit_steps(c, dt) + pad;
  for(int i = 1; i < 4; ++i) dims[i] = 2 * dimensions(i - 1);

  return dims;
}

std::vector<const_DotRange> AIM::Grid::box_contents_map(
    const DotVector &dots) const
{
  std::vector<const_DotRange> boxes(num_gridpoints);
  for(int box_idx = 0; box_idx < static_cast<int>(boxes.size()); ++box_idx) {
    auto NearestGridpoint = [=](const QuantumDot &qd) {
      return associated_grid_index(qd.position()) == box_idx;
    };

    auto begin = std::find_if(dots.begin(), dots.end(), NearestGridpoint);
    auto end = std::find_if_not(begin, dots.end(), NearestGridpoint);
    boxes.at(box_idx) = std::make_pair(begin, end);
  }

  return boxes;
}

std::vector<int> AIM::Grid::expansion_indices(const int grid_index) const
{
  Eigen::Vector3i origin = idx_to_coord(grid_index);
  std::vector<int> indices(std::pow(order_ + 1, 3));

  int idx = 0;
  for(int nx = 0; nx <= order_; ++nx) {
    for(int ny = 0; ny <= order_; ++ny) {
      for(int nz = 0; nz <= order_; ++nz) {
        const Eigen::Vector3i delta(Math::grid_sequence(nx),
                                    Math::grid_sequence(ny),
                                    Math::grid_sequence(nz));
        const int grid_idx = coord_to_idx(origin + delta);

        indices.at(idx++) = grid_idx;
      }
    }
  }

  return indices;
}

std::vector<AIM::Grid::ipair_t> AIM::Grid::nearfield_pairs(
    const int border, const DotVector &dots) const
{
  std::vector<ipair_t> nf;
  auto mapping = box_contents_map(dots);

  const int bound = order_ + border;

  for(int src_idx = 0; src_idx < num_gridpoints; ++src_idx) {
    Eigen::Vector3i s = idx_to_coord(src_idx);
    if(mapping[src_idx].first == mapping[src_idx].second) continue;

    for(int x = std::max(0, s(0) - bound);
        x <= std::min(dimensions(0) - 1, s(0) + bound); ++x) {
      for(int y = std::max(0, s(1) - bound);
          y <= std::min(dimensions(1) - 1, s(1) + bound); ++y) {
        for(int z = std::max(0, s(2) - bound);
            z <= std::min(dimensions(2) - 1, s(2) + bound); ++z) {
          const auto obs_idx = coord_to_idx({x, y, z});
          if(obs_idx >= src_idx &&
             mapping[obs_idx].first != mapping[obs_idx].second) {
            nf.emplace_back(src_idx, obs_idx);
          }
        }
      }
    }
  }

  nf.shrink_to_fit();
  return nf;
}

void AIM::Grid::sort_points_on_boxidx(DotVector &dots) const
{
  auto grid_comparitor = [&](const QuantumDot &q1, const QuantumDot &q2) {
    return associated_grid_index(q1.position()) <
           associated_grid_index(q2.position());
  };

  std::stable_sort(dots.begin(), dots.end(), grid_comparitor);
}

boost::multi_array<int, 2> AIM::Grid::expansion_index_table() const
{
  boost::multi_array<int, 2> neighborhood(
      boost::extents[num_gridpoints][std::pow(order_ + 1, 3)]);

  for(int i = 0; i < num_gridpoints; ++i) {
    auto indices = expansion_indices(i);
    std::copy(indices.begin(), indices.end(), &neighborhood[i][0]);
  }

  return neighborhood;
}
