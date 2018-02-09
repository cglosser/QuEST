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
  // Forget about last dimension; not needed in this scheme
  std::array<int, 5> extents = {table_dimensions[0], table_dimensions[1], 2,
                                table_dimensions[2], 3};
  source_table.resize(extents);
  obs_table.resize(extents);
  propagation_table = make_propagation_table();
}

void AIM::Nearfield::fill_source_table(const int step)
{
  using namespace Expansions::enums;

  const int wrapped_step = step % table_dimensions[0];

  for(int t = 0; t < table_dimensions[0]; ++t) {
    for(size_t p = 0; p < neighbors.size(); ++p) {
      size_t box1, box2;
      std::tie(box1, box2) = neighbors[p];

      for(size_t e = 0; e < expansion_table.shape()[1]; ++e) {
        Eigen::Map<Eigen::Vector3cd> grid_field1(
            &source_table[wrapped_step][p][0][e][0]);

        DotVector::const_iterator start, end;
        std::tie(start, end) = mapping[box1];
        for(auto d = start; d != end; ++d) {
          const auto dot_idx = std::distance(dots->begin(), d);
          grid_field1 += expansion_table[dot_idx][e].d0 *
                         (*dots)[dot_idx].dipole() *
                         history->array[dot_idx][step][0][RHO_01];
        }

        Eigen::Map<Eigen::Vector3cd> grid_field2(
            &source_table[wrapped_step][p][1][e][0]);

        std::tie(start, end) = mapping[box2];
        for(auto d = start; d != end; ++d) {
          const auto dot_idx = std::distance(dots->begin(), d);
          grid_field1 += expansion_table[dot_idx][e].d0 *
                         (*dots)[dot_idx].dipole() *
                         history->array[dot_idx][step][0][RHO_01];
        }
      }
    }
  }
}

void AIM::Nearfield::propagate(const int step)
{
  using Mat3d = Eigen::Matrix<cmplx, Eigen::Dynamic, 3>;
  const auto wrapped_step = step % table_dimensions[0];

  for(int t = 1; t < table_dimensions[0]; ++t) {
    const int wrap = std::max(step - t, 0) % table_dimensions[0];
    for(size_t p = 0; p < neighbors.size(); ++p) {
      Eigen::Map<Eigen::MatrixXcd> mat(&propagation_table[t][p][0][0],
                                       table_dimensions[2],
                                       table_dimensions[2]);

      Eigen::Map<Mat3d> src1(&source_table[wrap][p][0][0][0],
                             table_dimensions[2], 3);
      Eigen::Map<Mat3d> obs1(&obs_table[wrapped_step][p][0][0][0],
                             table_dimensions[2], 3);

      obs1 += mat * src1;

      if(neighbors[p].first == neighbors[p].second) continue;

      Eigen::Map<Mat3d> src2(&source_table[wrap][p][0][0][0],
                             table_dimensions[2], 3);
      Eigen::Map<Mat3d> obs2(&obs_table[wrapped_step][p][0][0][0],
                             table_dimensions[2], 3);

      obs2 += mat * src2;
    }
  }
}

void AIM::Nearfield::fill_results_table(const int step)
{
  using namespace Expansions::enums;

  const int wrapped_step = step % table_dimensions[0];
  results = 0;

  for(size_t p = 0; p < neighbors.size(); ++p) {
    size_t box1, box2;
    std::tie(box1, box2) = neighbors[p];

    for(size_t e = 0; e < expansion_table.shape()[1]; ++e) {
      Eigen::Map<Eigen::Vector3cd> grid_field1(
          &source_table[wrapped_step][p][0][e][0]);

      DotVector::const_iterator start, end;
      std::tie(start, end) = mapping[box1];
      for(auto d = start; d != end; ++d) {
        const auto dot_idx = std::distance(dots->begin(), d);
        results(dot_idx) += expansion_table[dot_idx][e].d0 *
                            grid_field1.dot((*dots)[dot_idx].dipole());
      }

      Eigen::Map<Eigen::Vector3cd> grid_field2(
          &source_table[wrapped_step][p][1][e][0]);

      std::tie(start, end) = mapping[box2];
      for(auto d = start; d != end; ++d) {
        const auto dot_idx = std::distance(dots->begin(), d);
        results(dot_idx) += expansion_table[dot_idx][e].d0 *
                            grid_field2.dot((*dots)[dot_idx].dipole());
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

        for(int t = 1; t < table_dimensions[0]; ++t) {
          const auto polynomial_idx = static_cast<int>(ceil(t - arg));
          if(0 <= polynomial_idx && polynomial_idx <= interp_order) {
            interp.evaluate_derivative_table_at_x(split_arg.second, dt);
            prop[t][p][e2][e1] =
                interp.evaluations[0][polynomial_idx] / normalization(dr);
          }
        }
      }
    }
  }

  return prop;
}
