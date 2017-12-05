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
                const Eigen::Vector3i &shift)
    : dimensions(dimensions),
      spacing(spacing),
      dots(nullptr),
      expansion_order(0)
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

Eigen::Vector3i AIM::Grid::grid_coordinate(const Eigen::Vector3d &coord) const
{
  return (coord.array() / spacing).cast<int>();
}

size_t AIM::Grid::coord_to_idx(const Eigen::Vector3i &coord) const
{
  return coord(2) + dimensions(2) * (coord(1) + dimensions(1) * coord(0));
}

Eigen::Vector3i AIM::Grid::idx_to_coord(size_t idx) const
{
  const int nynz = dimensions(1) * dimensions(2);
  const int x = idx / nynz;
  idx -= x * nynz;
  const int y = idx / dimensions(2);
  const int z = idx % dimensions(2);

  return Eigen::Vector3i(x, y, z);
}

Eigen::Vector3d AIM::Grid::spatial_coord_of_box(const size_t box_id) const
{
  return (idx_to_coord(box_id) + bounds.col(0).matrix())
             .array()
             .cast<double>() *
         spacing;
}

std::vector<size_t> AIM::Grid::expansion_box_indices(const Eigen::Vector3d &pos,
                                                     const int order) const
{
  Eigen::Vector3i origin = grid_coordinate(pos);
  std::vector<size_t> indices(std::pow(order + 1, 3));

  size_t idx = 0;
  for(int nx = 0; nx <= order; ++nx) {
    for(int ny = 0; ny <= order; ++ny) {
      for(int nz = 0; nz <= order; ++nz) {
        const Eigen::Vector3i delta(grid_sequence(nx), grid_sequence(ny),
                                    grid_sequence(nz));
        const size_t grid_idx = coord_to_idx(origin + delta);

        indices.at(idx++) = grid_idx;
      }
    }
  }

  return indices;
}

void AIM::Grid::sort_points_on_boxidx() const
{
  auto grid_comparitor = [&](const QuantumDot &q1, const QuantumDot &q2) {
    return coord_to_idx(grid_coordinate(q1.position())) <
           coord_to_idx(grid_coordinate(q2.position()));
  };

  std::stable_sort(dots->begin(), dots->end(), grid_comparitor);
}
