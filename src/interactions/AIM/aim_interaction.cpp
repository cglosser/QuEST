#include "aim_interaction.h"

AIM::AimInteraction::AimInteraction(const int interp_order,
                                    const Grid &grid,
                                    normalization::SpatialNorm normalization)
    : AimInteraction(nullptr,
                     nullptr,
                     interp_order,
                     1,
                     1,
                     grid,
                     Expansions::ExpansionTable(),
                     nullptr,
                     normalization)
{
}

AIM::AimInteraction::AimInteraction(
    const std::shared_ptr<const DotVector> dots,
    const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history,
    const int interp_order,
    const double c0,
    const double dt,
    const Grid &grid,
    const Expansions::ExpansionTable &expansion_table,
    Expansions::ExpansionFunction expansion_function,
    normalization::SpatialNorm normalization)
    : HistoryInteraction(dots, history, interp_order, c0, dt),
      grid(grid),
      expansion_table(expansion_table),
      expansion_function(std::move(expansion_function)),
      normalization(std::move(normalization)),

      // FFT stuff
      circulant_dimensions(grid.circulant_shape(c0, dt, interp_order)),
      fourier_table(circulant_fourier_table()),
      source_table(spacetime::make_vector3d<cmplx>(circulant_dimensions)),
      obs_table(spacetime::make_vector3d<cmplx>(circulant_dimensions)),
      spatial_vector_transforms(spatial_fft_plans()),

      // Nearfield stuff
      nf_pairs(grid.nearfield_pairs(2)),
      nf_matrices(boost::extents[nf_pairs.size()][circulant_dimensions[0]]
                                [expansion_table.shape()[1]]
                                [expansion_table.shape()[1]]),
      nf_source_table(expansion_table.shape()[1], 3),
      nf_obs_table(expansion_table.shape()[1], 3),
      nf_correction(spacetime::make_vector3d<cmplx>(circulant_dimensions))
{
  std::fill(source_table.data(),
            source_table.data() + source_table.num_elements(), cmplx(0, 0));
  std::fill(obs_table.data(), obs_table.data() + obs_table.num_elements(),
            cmplx(0, 0));

  fill_nearfield_matrices();
}

const Interaction::ResultArray &AIM::AimInteraction::evaluate(const int step)
{
  fill_source_table(step);
  propagate(step);
  fill_results_table(step);
  return results;
}

void AIM::AimInteraction::fill_source_table(const int step)
{
  using namespace Expansions::enums;

  const int wrapped_step = step % circulant_dimensions[0];
  const auto p = &source_table[wrapped_step][0][0][0][0];
  std::fill(p, p + 3 * 8 * grid.num_gridpoints, cmplx(0, 0));

  for(auto dot_idx = 0u; dot_idx < expansion_table.shape()[0]; ++dot_idx) {
    for(auto expansion_idx = 0u; expansion_idx < expansion_table.shape()[2];
        ++expansion_idx) {
      const Expansions::Expansion &e = expansion_table[dot_idx][expansion_idx];
      Eigen::Vector3i coord = grid.idx_to_coord(e.index);

      Eigen::Map<Eigen::Vector3cd> grid_field(
          &source_table[wrapped_step][coord(0)][coord(1)][coord(2)][0]);

      // This is the seam between what's stored in the History (density matrix
      // elements) and the electromagnetic source quantities. Ideally the AIM
      // code should not have knowledge of this to better encapsulate
      // "propagation," but this is good enough for now.
      Eigen::Vector3cd source_field = e.weights[D_0] *
                                      (*dots)[dot_idx].dipole() *
                                      history->array[dot_idx][step][0][RHO_01];
      grid_field += source_field;
    }
  }
}

void AIM::AimInteraction::propagate(const int step)
{
  const auto wrapped_step = step % circulant_dimensions[0];
  const auto nb = 8 * grid.num_gridpoints;
  const auto p = &source_table[wrapped_step][0][0][0][0];
  fftw_execute_dft(spatial_vector_transforms.forward,
                   reinterpret_cast<fftw_complex *>(p),
                   reinterpret_cast<fftw_complex *>(p));

  std::array<int, 5> front = {{wrapped_step, 0, 0, 0, 0}};
  Eigen::Map<Eigen::Array3Xcd> observers(&obs_table(front), 3, nb);
  observers = 0;

  for(int i = 1; i < circulant_dimensions[0]; ++i) {
    // If (step - i) runs "off the end", just propagate src[0][...]
    auto wrap = std::max(step - i, 0) % circulant_dimensions[0];

    Eigen::Map<Eigen::ArrayXcd> prop(&fourier_table[i][0][0][0], nb);
    Eigen::Map<Eigen::Array3Xcd> src(&source_table[wrap][0][0][0][0], 3, nb);

    // Use broadcasting to do the x, y, and z component propagation
    observers += src.rowwise() * prop.transpose();
  }

  fftw_execute_dft(spatial_vector_transforms.backward,
                   reinterpret_cast<fftw_complex *>(observers.data()),
                   reinterpret_cast<fftw_complex *>(observers.data()));
}

void AIM::AimInteraction::fill_results_table(const int step)
{
  results = 0;

  for(auto dot_idx = 0u; dot_idx < expansion_table.shape()[0]; ++dot_idx) {
    Eigen::Vector3cd total_field = Eigen::Vector3cd::Zero();
    for(auto expansion_idx = 0u; expansion_idx < expansion_table.shape()[1];
        ++expansion_idx) {
      const Expansions::Expansion &e = expansion_table[dot_idx][expansion_idx];
      Eigen::Vector3i coord = grid.idx_to_coord(e.index);
      total_field += expansion_function(
          obs_table, {{step, coord(0), coord(1), coord(2)}}, e);
      // Don't use a _wrapped_ step here; the expansion_function needs knowledge
      // of where it's being called in the complete timeline to accommodate
      // boundary conditions
    }
    results(dot_idx) += total_field.dot((*dots)[dot_idx].dipole());
  }
}

spacetime::vector<cmplx> AIM::AimInteraction::circulant_fourier_table()
{
  spacetime::vector<cmplx> g_mat(circulant_dimensions);

  const int num_gridpts = circulant_dimensions[1] * circulant_dimensions[2] *
                          circulant_dimensions[3];
  TransformPair circulant_plan = {
      fftw_plan_many_dft(3, &circulant_dimensions[1], circulant_dimensions[0],
                         reinterpret_cast<fftw_complex *>(g_mat.data()),
                         nullptr, 1, num_gridpts,
                         reinterpret_cast<fftw_complex *>(g_mat.data()),
                         nullptr, 1, num_gridpts, FFTW_FORWARD, FFTW_MEASURE),
      nullptr};

  fill_gmatrix_table(g_mat);

  // Transform the circulant vectors into their equivalently-diagonal
  // representation. Buckle up.

  fftw_execute(circulant_plan.forward);

  // This accounts for FFTW's *un*normalized transform -- it takes the least
  // amount of computational effort to put all of the normalizations here.

  Eigen::Map<Eigen::ArrayXcd> gs(g_mat.data(), g_mat.num_elements());
  gs /= num_gridpts;

  return g_mat;
}

void AIM::AimInteraction::fill_gmatrix_table(
    spacetime::vector<cmplx> &gmatrix_table) const
{  // Build the circulant vectors that define the G "matrices." Since the G
  // matrices are Toeplitz (and symmetric), they're uniquely determined by
  // their first row. The first row gets computed here then mirrored to make a
  // list of every circulant (and thus FFT-able) vector. This function needs to
  // accept a non-const reference to a spacetime::vector (instead of just
  // returning such an array) to play nice with FFTW and its workspaces.

  std::fill(gmatrix_table.data(),
            gmatrix_table.data() + gmatrix_table.num_elements(),
            cmplx(0.0, 0.0));

  Interpolation::UniformLagrangeSet interp(interp_order);

  for(int x = 0; x < grid.dimensions(0); ++x) {
    for(int y = 0; y < grid.dimensions(1); ++y) {
      for(int z = 0; z < grid.dimensions(2); ++z) {
        bool nearfield_point =
            (x <= interp_order) && (y <= interp_order) && (z <= interp_order);
        if(nearfield_point) continue;

        const size_t box_idx = grid.coord_to_idx(Eigen::Vector3i(x, y, z));
        const Eigen::Vector3d dr =
            grid.spatial_coord_of_box(box_idx) - grid.spatial_coord_of_box(0);

        const double arg = dr.norm() / (c0 * dt);
        const std::pair<int, double> split_arg = split_double(arg);

        for(int t = 1; t < circulant_dimensions[0]; ++t) {
          const auto polynomial_idx = static_cast<int>(ceil(t - arg));
          if(0 <= polynomial_idx && polynomial_idx <= interp_order) {
            interp.evaluate_derivative_table_at_x(split_arg.second, dt);
            gmatrix_table[t][x][y][z] =
                interp.evaluations[0][polynomial_idx] / normalization(dr);
          }
        }
      }
    }
  }

  spacetime::fill_circulant_mirror(gmatrix_table);
}

TransformPair AIM::AimInteraction::spatial_fft_plans()
{
  // Set up FFTW plans to transform projected source distributions. Due to the
  // requirements of the circulant extension, these plans perform transforms of
  // length 2 n_{x,y,z} to accommodate the requisite zero padding. While they're
  // constructed to work on the head of `source_table` (that is, what would be
  // the I_0 source), the advanced FFTW interface allows them to stride forward
  // to equivalently transform the source currents at every timestep.

  constexpr int num_transforms = 3;
  constexpr int transform_rank = 3;
  constexpr int dist_between_elements = 3;
  constexpr int dist_between_transforms = 1;

  auto make_plan = [&](const int sign) {
    return fftw_plan_many_dft(
        transform_rank, &circulant_dimensions[1], num_transforms,
        reinterpret_cast<fftw_complex *>(source_table.data()), nullptr,
        dist_between_elements, dist_between_transforms,
        reinterpret_cast<fftw_complex *>(source_table.data()), nullptr,
        dist_between_elements, dist_between_transforms, sign, FFTW_MEASURE);
  };

  return {make_plan(FFTW_FORWARD), make_plan(FFTW_BACKWARD)};
}

void AIM::AimInteraction::fill_nearfield_matrices()
{
  std::fill(nf_matrices.data(), nf_matrices.data() + nf_matrices.num_elements(),
            0.0);
  Interpolation::UniformLagrangeSet interp(interp_order);

  for(auto p = 0u; p < nf_pairs.size(); ++p) {
    const auto src_indices = grid.expansion_indices(nf_pairs[p].first);
    const auto obs_indices = grid.expansion_indices(nf_pairs[p].second);

    for(auto i = 0u; i < src_indices.size(); ++i) {
      for(auto j = 0u; j < obs_indices.size(); ++j) {
        Eigen::Vector3d dr = grid.spatial_coord_of_box(src_indices[i]) -
                             grid.spatial_coord_of_box(obs_indices[j]);
        const double arg = dr.norm() / (c0 * dt);
        const std::pair<int, double> split_arg = split_double(arg);
        for(int t = 1; t < circulant_dimensions[0]; ++t) {
          const auto polynomial_idx = static_cast<int>(ceil(t - arg));
          if(0 <= polynomial_idx && polynomial_idx <= interp_order) {
            interp.evaluate_derivative_table_at_x(split_arg.second, dt);
            nf_matrices[p][t][i][j] =
                interp.evaluations[0][polynomial_idx] / normalization(dr);
          }
        }
      }
    }
  }
}

void AIM::AimInteraction::evaluate_nearfield(const int step)
{
  const int wrapped_step = step % circulant_dimensions[0];
  std::fill(nf_correction.data(),
            nf_correction.data() + nf_correction.num_elements(), cmplx(0, 0));

  for(auto p = 0u; p < nf_pairs.size(); ++p) {
    const auto src_indices = grid.expansion_indices(nf_pairs[p].first);
    const auto obs_indices = grid.expansion_indices(nf_pairs[p].second);

    nf_obs_table.setZero();

    for(int t = 1; t < circulant_dimensions[0]; ++t) {
      const auto wrap = std::max(step - t, 0) % circulant_dimensions[0];

      // Grab the grid values for the source box and condense them into a "dense
      // vector"...
      for(auto s = 0u; s < src_indices.size(); ++s) {
        Eigen::Array3i coord = grid.idx_to_coord(src_indices[s]);
        nf_source_table.row(s) = Eigen::Map<Eigen::Vector3cd>(
            &source_table[wrap][coord[0]][coord[1]][coord[2]][0]);
      }

      Eigen::Map<Eigen::MatrixXcd> g(&nf_matrices[p][t][0][0],
                                     nf_matrices.shape()[2],
                                     nf_matrices.shape()[3]);

      // ...propagate them with a small, dense propagator...
      nf_obs_table.col(0) += g * nf_source_table.col(0);
      nf_obs_table.col(1) += g * nf_source_table.col(1);
      nf_obs_table.col(2) += g * nf_source_table.col(2);
    }

    //...and put them back in a global-sized table.
    for(auto u = 0u; u < obs_indices.size(); ++u) {
      Eigen::Array3i coord = grid.idx_to_coord(obs_indices[u]);
      Eigen::Map<Eigen::Vector3cd> dest(
          &nf_correction[wrapped_step][coord[0]][coord[1]][coord[2]][0]);
      dest = nf_obs_table.row(u);
    }

    // Same box; need to avoid double-counting
    if(nf_pairs[p].first == nf_pairs[p].second) continue;

    // Now do the same as above; just swap src and obs
    nf_obs_table.setZero();

    for(int t = 1; t < circulant_dimensions[0]; ++t) {
      const auto wrap = std::max(step - t, 0) % circulant_dimensions[0];

      for(auto s = 0u; s < obs_indices.size(); ++s) {
        Eigen::Array3i coord = grid.idx_to_coord(obs_indices[s]);
        nf_source_table.row(s) = Eigen::Map<Eigen::Vector3cd>(
            &source_table[wrap][coord[0]][coord[1]][coord[2]][0]);
      }

      Eigen::Map<Eigen::MatrixXcd> g(&nf_matrices[p][t][0][0],
                                     nf_matrices.shape()[2],
                                     nf_matrices.shape()[3]);

      nf_obs_table.col(0) += g * nf_source_table.col(0);
      nf_obs_table.col(1) += g * nf_source_table.col(1);
      nf_obs_table.col(2) += g * nf_source_table.col(2);
    }

    for(auto u = 0u; u < obs_indices.size(); ++u) {
      Eigen::Array3i coord = grid.idx_to_coord(src_indices[u]);
      Eigen::Map<Eigen::Vector3cd> dest(
          &nf_correction[wrapped_step][coord[0]][coord[1]][coord[2]][0]);
      dest = nf_obs_table.row(u);
    }
  }
}
