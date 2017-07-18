#include "aim_interaction.h"

AIM::Grid::Grid(const Eigen::Array3d &spacing,
                const std::shared_ptr<DotVector> &dots)
    : spacing(spacing), dots(dots), bounds(calculate_bounds())
{
  num_boxes = bounds.col(1) - bounds.col(0) + 1;
  max_diagonal = (num_boxes.cast<double>() * spacing).matrix().norm();
  boxes.resize(num_boxes.prod());

  sort_points_on_boxidx();
  map_points_to_boxes();
}

AIM::Grid::BoundsArray AIM::Grid::calculate_bounds() const
{
  BoundsArray b;
  b.setZero();

  for(const auto &qdot : *dots) {
    Eigen::Vector3i grid_coord = grid_coordinate(qdot.position());

    b.col(0) = grid_coord.array().min(b.col(0));
    b.col(1) = grid_coord.array().max(b.col(1));
  }

  return b;
}

Eigen::Vector3i AIM::Grid::grid_coordinate(const Eigen::Vector3d &coord) const
{
  return floor(coord.cwiseQuotient(spacing.matrix()).array()).cast<int>();
}

size_t AIM::Grid::coord_to_idx(const Eigen::Vector3i &coord) const
{
  Eigen::Vector3i shifted(coord - bounds.col(0).matrix());

  return shifted(0) + num_boxes(0) * (shifted(1) + num_boxes(1) * shifted(2));
}

Eigen::Vector3i AIM::Grid::idx_to_coord(size_t idx) const
{
  const int nxny = num_boxes(0) * num_boxes(1);
  const int z = idx / nxny;
  idx -= z * nxny;
  const int y = idx / num_boxes(0);
  const int x = idx % num_boxes(0);

  return Eigen::Vector3i(x, y, z);
}

void AIM::Grid::sort_points_on_boxidx() const
{
  auto grid_comparitor = [&](const QuantumDot &q1, const QuantumDot &q2) {
    return coord_to_idx(grid_coordinate(q1.position())) <
           coord_to_idx(grid_coordinate(q2.position()));
  };

  std::stable_sort(dots->begin(), dots->end(), grid_comparitor);
}

void AIM::Grid::map_points_to_boxes()
{
  for(size_t box_idx = 0; box_idx < boxes.size(); ++box_idx) {
    auto IsInBox = [=](const QuantumDot &qd) {
      return coord_to_idx(grid_coordinate(qd.position())) == box_idx;
    };

    auto begin = std::find_if(dots->begin(), dots->end(), IsInBox);
    auto end = std::find_if_not(begin, dots->end(), IsInBox);
    boxes.at(box_idx) = std::make_pair(begin, end);
  }
}

AIM::AimInteraction::AimInteraction(const std::shared_ptr<DotVector> &dots,
                                    const Eigen::Vector3d &spacing)
    : Interaction(dots), grid(spacing, dots)
{
}
