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

std::vector<AIM::Grid::ipair_t> AIM::Grid::nearfield_pairs(
    const int thresh) const
{
  std::vector<ipair_t> nf;

  const int bound = expansion_order + thresh;

  for(auto src_idx = 0u; src_idx < num_gridpoints; ++src_idx) {
    Eigen::Vector3i s = idx_to_coord(src_idx);
    for(int x = std::max(0, s(0) - bound);
        x <= std::min(dimensions(0) - 1, s(0) + bound); ++x) {
      for(int y = std::max(0, s(1) - bound);
          y <= std::min(dimensions(1) - 1, s(1) + bound); ++y) {
        for(int z = std::max(0, s(2) - bound);
            z <= std::min(dimensions(2) - 1, s(2) + bound); ++z) {
          const auto obs_idx = coord_to_idx({x, y, z});
          if(obs_idx >= src_idx) {
            nf.emplace_back(src_idx, obs_idx);
          }
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
