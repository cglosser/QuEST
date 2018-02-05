#include "aim_interaction.h"

AIM::AimInteraction::AimInteraction(const int interp_order,
                                    const Grid grid,
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
    const Grid grid,
    const Expansions::ExpansionTable &expansion_table,
    Expansions::ExpansionFunction expansion_function,
    normalization::SpatialNorm normalization)
    : HistoryInteraction(dots, history, interp_order, c0, dt),
      grid(std::move(grid)),
      expansion_table(expansion_table),
      expansion_function(std::move(expansion_function)),
      normalization(std::move(normalization)),

      // FFT stuff
      toeplitz_dimensions(grid.toeplitz_shape(c0, dt, interp_order)),
      circulant_dimensions(grid.circulant_shape(c0, dt, interp_order)),
      fourier_table(circulant_fourier_table()),
      source_table(spacetime::make_vector3d<cmplx>(circulant_dimensions)),
      source_table_fft(spacetime::make_vector3d<cmplx>(circulant_dimensions)),
      obs_table(spacetime::make_vector3d<cmplx>(circulant_dimensions)),
      obs_table_fft(spacetime::make_vector3d<cmplx>(circulant_dimensions)),
      spatial_vector_transforms(spatial_fft_plans()),

      // Nearfield stuff
      nf_pairs(grid.full_nearfield_pairs(1)),
      nf_mats_sparse(generate_nearfield_mats()),
      nf_correction(grid.num_gridpoints, 3),
      nf_workspace(grid.num_gridpoints, 3)
{
  auto clear = [](auto &table) {
    std::fill(table.data(), table.data() + table.num_elements(), cmplx(0, 0));
  };

  clear(source_table);
  clear(source_table_fft);
  clear(obs_table);
  clear(obs_table_fft);
}

const Interaction::ResultArray &AIM::AimInteraction::evaluate(const int step)
{
  fill_source_table(step);
  propagate(step);
  fill_results_table(step);

  evaluate_nearfield(step);

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
      Eigen::Vector3cd source_field = e.d0 * (*dots)[dot_idx].dipole() *
                                      history->array[dot_idx][step][0][RHO_01];
      grid_field += source_field;
    }
  }
}

void AIM::AimInteraction::propagate(const int step)
{
  const auto wrapped_step = step % circulant_dimensions[0];
  const auto nb = 8 * grid.num_gridpoints;
  const std::array<int, 5> front = {{wrapped_step, 0, 0, 0, 0}};

  const auto s_ptr = &source_table(front), s_ptr_fft = &source_table_fft(front);
  fftw_execute_dft(spatial_vector_transforms.forward,
                   reinterpret_cast<fftw_complex *>(s_ptr),
                   reinterpret_cast<fftw_complex *>(s_ptr_fft));

  Eigen::Map<Eigen::Array3Xcd> observers(&obs_table_fft(front), 3, nb);
  observers = 0;

  for(int i = 1; i < circulant_dimensions[0]; ++i) {
    // If (step - i) runs "off the end", just propagate src[0][...]
    auto wrap = std::max(step - i, 0) % circulant_dimensions[0];

    Eigen::Map<Eigen::ArrayXcd> prop(&fourier_table[i][0][0][0], nb);
    Eigen::Map<Eigen::Array3Xcd> src(&source_table_fft[wrap][0][0][0][0], 3,
                                     nb);

    // Use broadcasting to do the x, y, and z component propagation
    observers += src.rowwise() * prop.transpose();
  }

  const auto o_ptr = &obs_table(front), o_ptr_fft = &obs_table_fft(front);
  fftw_execute_dft(spatial_vector_transforms.backward,
                   reinterpret_cast<fftw_complex *>(o_ptr_fft),
                   reinterpret_cast<fftw_complex *>(o_ptr));
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
        const size_t box_idx = grid.coord_to_idx({x, y, z});

        if(box_idx == 0) continue;

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

void AIM::AimInteraction::evaluate_nearfield(const int step)
{
  nf_correction.setZero();

  for(int t = 1; t < circulant_dimensions[0]; ++t) {
    const auto wrap = std::max(step - t, 0) % circulant_dimensions[0];

    int i = 0;
    for(int x = 0; x < toeplitz_dimensions[1]; ++x) {
      for(int y = 0; y < toeplitz_dimensions[2]; ++y) {
        for(int z = 0; z < toeplitz_dimensions[3]; ++z) {
          nf_workspace.row(i++) =
              Eigen::Map<Eigen::Vector3cd>(&source_table[wrap][x][y][z][0]);
          // This is necessary to "dodge" the circulant holes of source_table
          // so as to simplify the (sparse) matrix multiplication
        }
      }
    }

    nf_correction += nf_mats_sparse[t] * nf_workspace;
  }

  const int wrapped_step = step % circulant_dimensions[0];
  int i = 0;
  for(int x = 0; x < toeplitz_dimensions[1]; ++x) {
    for(int y = 0; y < toeplitz_dimensions[2]; ++y) {
      for(int z = 0; z < toeplitz_dimensions[3]; ++z) {
        Eigen::Map<Eigen::Vector3cd>(&obs_table[wrapped_step][x][y][z][0]) -=
            nf_correction.row(i++);
      }
    }
  }
}

std::vector<Eigen::SparseMatrix<cmplx>>
AIM::AimInteraction::generate_nearfield_mats()
{
  using Triplet = Eigen::Triplet<cmplx>;

  Interpolation::UniformLagrangeSet interp(interp_order);

  std::vector<std::vector<Triplet>> coefficients(circulant_dimensions[0]);
  for(const auto &p : nf_pairs) {
    const Eigen::Vector3d dr = grid.spatial_coord_of_box(p.first) -
                               grid.spatial_coord_of_box(p.second);

    const double arg = dr.norm() / (c0 * dt);
    const auto split_arg = split_double(arg);

    for(int t = 0; t < circulant_dimensions[0]; ++t) {
      const auto polynomial_idx = static_cast<int>(ceil(t - arg));
      if(0 <= polynomial_idx && polynomial_idx <= interp_order) {
        interp.evaluate_derivative_table_at_x(split_arg.second, dt);
        const auto val =
            interp.evaluations[0][polynomial_idx] / normalization(dr);
        coefficients[t].push_back(Triplet(p.first, p.second, val));
      }
    }
  }

  std::vector<Eigen::SparseMatrix<cmplx>> mats(
      circulant_dimensions[0],
      Eigen::SparseMatrix<cmplx>(grid.num_gridpoints, grid.num_gridpoints));

  for(int i = 0; i < circulant_dimensions[0]; ++i) {
    mats.at(i).setFromTriplets(coefficients.at(i).begin(),
                               coefficients.at(i).end());
  }

  return mats;
}
