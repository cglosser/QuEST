#include "nearfield.h"

AIM::Nearfield::Nearfield(
    const std::shared_ptr<const DotVector> dots,
    const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history,
    const int interp_order,
    const int border,
    const double c0,
    const double dt,
    const Grid grid,
    const Expansions::ExpansionTable &expansion_table,
    Expansions::ExpansionFunction expansion_function,
    normalization::SpatialNorm normalization)
    : AimBase(dots,
              history,
              interp_order,
              c0,
              dt,
              grid,
              expansion_table,
              expansion_function,
              normalization,
              grid.nearfield_shape(c0, dt, interp_order, border)),
      mapping(grid.box_contents_map()),
      neighbors(grid.nearfield_pairs(border))
{
  // The order here is time, pair, point set, expansion index, vector dimension
  // in accordance with the FIELD_AXIS_LABEL enum. Each matrix product maps a
  // set of E points into a set of E points -- these sets are what are indexed
  // by "point set" (0 corresponds to expansion points for the lower-indexed box
  // and 1 corresponds to the same for the equally or higher indexed box, as
  // determined by the structure of neighbors.
  // It's the biggest hack that the dimensions of these arrays are the same for
  // the NF and FF code. Sorry.
  field_table_dims = {table_dimensions[0], table_dimensions[1], 2,
                      table_dimensions[2], 3};
  source_table.resize(field_table_dims);
  obs_table.resize(field_table_dims);
  propagation_table = make_propagation_table();
}

void AIM::Nearfield::fill_source_table(const int step)
{
  using namespace Expansions::enums;

  const int wrapped_step = step % table_dimensions[0];
  spacetime::clear_time_slice(source_table, wrapped_step);

  for(int p = 0; p < field_table_dims[PAIRS]; ++p) {
    size_t box1, box2;
    std::tie(box1, box2) = neighbors[p];

    for(int e = 0; e < field_table_dims[EXPANSIONS]; ++e) {
      Eigen::Map<Eigen::Vector3cd> grid_field1(
          &source_table[wrapped_step][p][0][e][0]);

      DotVector::const_iterator start, end;
      std::tie(start, end) = mapping[box1];
      for(auto d = start; d != end; ++d) {
        const auto dot_idx = std::distance(dots->begin(), d);
        grid_field1 += expansion_table[dot_idx][e].d0 *
                       (*dots)[dot_idx].dipole() *
                       history->array_[dot_idx][step][0][RHO_01];
      }

      Eigen::Map<Eigen::Vector3cd> grid_field2(
          &source_table[wrapped_step][p][1][e][0]);

      std::tie(start, end) = mapping[box2];
      for(auto d = start; d != end; ++d) {
        const auto dot_idx = std::distance(dots->begin(), d);
        grid_field2 += expansion_table[dot_idx][e].d0 *
                       (*dots)[dot_idx].dipole() *
                       history->array_[dot_idx][step][0][RHO_01];
      }
    }
  }
}

void AIM::Nearfield::propagate(const int step)
{
  using Mat3d = Eigen::Matrix<cmplx, Eigen::Dynamic, 3, Eigen::RowMajor>;

  const auto wrapped_step = step % table_dimensions[0];
  spacetime::clear_time_slice(obs_table, wrapped_step);

  for(int t = 0; t < field_table_dims[STEPS]; ++t) {
    const int wrap = std::max(step - t, 0) % table_dimensions[0];
    for(int p = 0; p < field_table_dims[PAIRS]; ++p) {
      Eigen::Map<Eigen::MatrixXcd> mat(&propagation_table[t][p][0][0],
                                       table_dimensions[2],
                                       table_dimensions[2]);

      Eigen::Map<Mat3d> src1(&source_table[wrap][p][0][0][0],
                             table_dimensions[2], 3);
      Eigen::Map<Mat3d> obs1(&obs_table[wrapped_step][p][1][0][0],
                             table_dimensions[2], 3);

      obs1 += mat * src1;

      if(neighbors[p].first == neighbors[p].second) continue;

      Eigen::Map<Mat3d> src2(&source_table[wrap][p][1][0][0],
                             table_dimensions[2], 3);
      Eigen::Map<Mat3d> obs2(&obs_table[wrapped_step][p][0][0][0],
                             table_dimensions[2], 3);

      obs2 += mat.transpose() * src2;
    }
  }
}

void AIM::Nearfield::fill_results_table(const int step)
{
  using namespace Expansions::enums;
  results = 0;

  for(int p = 0; p < field_table_dims[PAIRS]; ++p) {
    size_t box1, box2;
    std::tie(box1, box2) = neighbors[p];

    for(int e = 0; e < field_table_dims[EXPANSIONS]; ++e) {
      DotVector::const_iterator start, end;
      std::tie(start, end) = mapping[box1];
      for(auto d = start; d != end; ++d) {
        const auto dot_idx = std::distance(dots->begin(), d);
        results(dot_idx) += (*dots)[dot_idx].dipole().dot(expansion_function(
            obs_table, {{step, p, 0, e}}, expansion_table[dot_idx][e]));
      }

      std::tie(start, end) = mapping[box2];
      for(auto d = start; d != end; ++d) {
        const auto dot_idx = std::distance(dots->begin(), d);
        results(dot_idx) += (*dots)[dot_idx].dipole().dot(expansion_function(
            obs_table, {{step, p, 1, e}}, expansion_table[dot_idx][e]));
      }
    }
  }
}

spacetime::vector<cmplx> AIM::Nearfield::make_propagation_table() const
{
  spacetime::vector<cmplx> prop(table_dimensions);
  std::fill(prop.data(), prop.data() + prop.num_elements(), cmplx(0.0, 0.0));

  Interpolation::UniformLagrangeSet interp(interp_order);

  for(size_t p = 0; p < neighbors.size(); ++p) {
    const auto expansion1 = grid.expansion_indices(neighbors[p].first);
    const auto expansion2 = grid.expansion_indices(neighbors[p].second);

    for(size_t e1 = 0; e1 < expansion_table.shape()[1]; ++e1) {
      for(size_t e2 = 0; e2 < expansion_table.shape()[1]; ++e2) {
        if(expansion1[e1] == expansion2[e2]) continue;
        Eigen::Vector3d dr = grid.spatial_coord_of_box(expansion1[e1]) -
                             grid.spatial_coord_of_box(expansion2[e2]);

        const double arg = dr.norm() / (c0 * dt);
        const std::pair<int, double> split_arg = split_double(arg);

        for(int t = 0; t < table_dimensions[0]; ++t) {
          const auto polynomial_idx = static_cast<int>(ceil(t - arg));
          if(0 <= polynomial_idx && polynomial_idx <= interp_order) {
            interp.evaluate_derivative_table_at_x(split_arg.second, dt);
            prop[t][p][e1][e2] =
                interp.evaluations[0][polynomial_idx] / normalization(dr);
          }
        }
      }
    }
  }

  return prop;
}
