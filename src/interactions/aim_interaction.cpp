#include "aim_interaction.h"

AIM::Grid::Grid(const Eigen::Array3d &spacing,
                const std::shared_ptr<DotVector> &dots)
    : spacing(spacing), dots(dots), bounds(calculate_bounds())
{
  dimensions = bounds.col(1) - bounds.col(0) + 1;
  num_boxes = dimensions.prod();
  max_diagonal = (dimensions.cast<double>() * spacing).matrix().norm();
  boxes.resize(dimensions.prod());

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

  return shifted(0) + dimensions(0) * (shifted(1) + dimensions(1) * shifted(2));
}

Eigen::Vector3i AIM::Grid::idx_to_coord(size_t idx) const
{
  const int nxny = dimensions(0) * dimensions(1);
  const int z = idx / nxny;
  idx -= z * nxny;
  const int y = idx / dimensions(0);
  const int x = idx % dimensions(0);

  return Eigen::Vector3i(x, y, z);
}

Eigen::Vector3d AIM::Grid::spatial_coord_of_box(const size_t box_id) const
{
  const Eigen::Vector3d r =
      (idx_to_coord(box_id).cast<double>().array() * spacing);
  return r + bounds.col(0).cast<double>().matrix();
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
                                    const Eigen::Vector3d &spacing,
                                    const int interp_order,
                                    const double c,
                                    const double dt)
    : Interaction(dots),
      grid(spacing, dots),
      interp_order(interp_order),
      c(c),
      dt(dt)
{
}

Interaction::ResultArray &AIM::AimInteraction::evaluate(const int step)
{
  results = 0;
  return results;
}

std::vector<double> AIM::AimInteraction::g_matrix_row(const size_t step) const
{
  std::vector<double> row(grid.num_boxes, 0);

  Interpolation::UniformLagrangeSet interp(interp_order);
  for(size_t box_idx = 1; box_idx < grid.num_boxes; ++box_idx) {
    const Eigen::Vector3d dr =
        grid.spatial_coord_of_box(box_idx) - grid.spatial_coord_of_box(0);

    const double arg = dr.norm() / (c * dt);
    const std::pair<int, double> split_arg = split_double(arg);

    const int polynomial_idx = static_cast<int>(ceil(step - arg));

    interp.evaluate_derivative_table_at_x(split_arg.second, dt);

    if(0 <= polynomial_idx && polynomial_idx <= interp_order) {
      row.at(box_idx) = interp.evaluations[0][polynomial_idx];
    }
  }

  return row;
}
