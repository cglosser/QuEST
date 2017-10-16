#include "aim_interaction.h"

AIM::AimInteraction::AimInteraction(
    const std::shared_ptr<const DotVector> &dots,
    const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> &history,
    const std::shared_ptr<Propagation::RotatingFramePropagator> &propagator,
    const int interp_order,
    const double c0,
    const double dt,
    const Grid &grid,
    const Array<Expansion> &expansion_table,
    const normalization::SpatialNorm &normalization)
    : HistoryInteraction(dots, history, propagator, interp_order, c0, dt),
      grid(grid),
      expansion_table(expansion_table),
      normalization(normalization),
      max_transit_steps(grid.max_transit_steps(c0, dt)),
      circulant_dimensions(grid.circulant_shape(c0, dt)),
      fourier_table(circulant_fourier_table()),
      source_table(circulant_dimensions),
      obs_table(circulant_dimensions),
      spatial_transforms(spatial_fft_plans())
{
  std::fill(source_table.data(),
            source_table.data() + source_table.num_elements(), cmplx(0, 0));
  std::fill(obs_table.data(), obs_table.data() + obs_table.num_elements(),
            cmplx(0, 0));
}

const Interaction::ResultArray &AIM::AimInteraction::evaluate(const int step)
{
  const int nb = 8 * grid.num_boxes;
  Eigen::Map<Eigen::ArrayXcd> obs_vec(&obs_table[step][0][0][0], nb);
  obs_vec = 0;

  for(int i = 0; i < step; ++i) {
    Eigen::Map<Eigen::ArrayXcd> prop(&fourier_table[step - i][0][0][0], nb);
    Eigen::Map<Eigen::ArrayXcd> src(&source_table[i][0][0][0], nb);

    obs_vec += prop * src;
  }

  fftw_execute_dft(spatial_transforms.backward,
                   reinterpret_cast<fftw_complex *>(&obs_table[step][0][0][0]),
                   reinterpret_cast<fftw_complex *>(&obs_table[step][0][0][0]));

  results = ResultArray(grid.num_boxes);

  int i = 0;
  for(auto x = 0l; x < grid.dimensions(0); ++x) {
    for(auto y = 0l; y < grid.dimensions(1); ++y) {
      for(auto z = 0l; z < grid.dimensions(2); ++z) {
        results(i++) = obs_table[step][x][y][z];
      }
    }
  }

  return results;
}

void AIM::AimInteraction::fill_source_table(const int step)
{
  for(auto dot_idx = 0u; dot_idx < expansion_table.shape()[0]; ++dot_idx) {
    for(auto expansion_idx = 0u; expansion_idx < expansion_table.shape()[1];
        ++expansion_idx) {
      const Expansion &e = expansion_table[dot_idx][expansion_idx];
      Eigen::Vector3i coord = grid.idx_to_coord(e.index);

      // This is the seam between what's stored in the History (density matrix
      // elements) and the electromagnetic source quantities. Ideally the AIM
      // code should not have knowledge of this to better encapsulate
      // "propagation," but this is good enough for now.
      source_table[step][coord(0)][coord(1)][coord(2)] =
          e.weight * history->array[dot_idx][step][0][1];
    }
  }
}

SpacetimeVector<cmplx> AIM::AimInteraction::circulant_fourier_table()
{
  SpacetimeVector<cmplx> g_mat(circulant_dimensions);

  const int num_gridpts = circulant_dimensions[1] * circulant_dimensions[2] *
                          circulant_dimensions[3];
  TransformPair circulant_plan = {
      fftw_plan_many_dft(3, &circulant_dimensions[1], circulant_dimensions[0],
                         reinterpret_cast<fftw_complex *>(g_mat.data()),
                         nullptr, 1, num_gridpts,
                         reinterpret_cast<fftw_complex *>(g_mat.data()),
                         nullptr, 1, num_gridpts, FFTW_FORWARD, FFTW_MEASURE),
      nullptr};

  std::fill(g_mat.data(), g_mat.data() + g_mat.num_elements(), cmplx(0, 0));

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
    SpacetimeVector<cmplx> &gmatrix_table) const
{  // Build the circulant vectors that define the G "matrices." Since the G
  // matrices are Toeplitz (and symmetric), they're uniquely determined by
  // their first row. The first row gets computed here then mirrored to make a
  // list of every circulant (and thus FFT-able) vector. This function needs to
  // accept a non-const reference to a SpacetimeVector (instead of just
  // returning such an array) to play nice with FFTW and its workspaces.

  Interpolation::UniformLagrangeSet interp(interp_order);
  for(int x = 0; x < grid.dimensions(0); ++x) {
    for(int y = 0; y < grid.dimensions(1); ++y) {
      for(int z = 0; z < grid.dimensions(2); ++z) {
        const size_t box_idx = grid.coord_to_idx(Eigen::Vector3i(x, y, z));
        if(box_idx == 0) continue;

        const auto dr =
            grid.spatial_coord_of_box(box_idx) - grid.spatial_coord_of_box(0);

        const double arg = dr.norm() / (c0 * dt);
        const std::pair<int, double> split_arg = split_double(arg);

        for(int t = 1; t < circulant_dimensions[0]; ++t) {
          const int polynomial_idx = static_cast<int>(ceil(t - arg));
          if(0 <= polynomial_idx && polynomial_idx <= interp_order) {
            interp.evaluate_derivative_table_at_x(split_arg.second, dt);
            gmatrix_table[t][x][y][z] =
                interp.evaluations[0][polynomial_idx] / normalization(dr);
          }
        }
      }
    }
  }

  fill_circulant_mirror(gmatrix_table);
}

TransformPair AIM::AimInteraction::spatial_fft_plans()
{
  // Set up FFTW plans to transform projected source distributions. Due to the
  // requirements of the circulant extension, these plans perform transforms of
  // length 2 n_{x,y,z} to accommodate the requisite zero padding. While they're
  // constructed to work on the head of `source_table` (that is, what would be
  // the I_0 source), the advanced FFTW interface allows them to stride forward
  // to equivalently transform the source currents at every timestep.

  auto make_plan = [&](const int sign) {
    return fftw_plan_dft_3d(
        circulant_dimensions[1], circulant_dimensions[2],
        circulant_dimensions[3],
        reinterpret_cast<fftw_complex *>(source_table.data()),
        reinterpret_cast<fftw_complex *>(source_table.data()), sign,
        FFTW_MEASURE);
  };

  return {make_plan(FFTW_FORWARD), make_plan(FFTW_BACKWARD)};
}
