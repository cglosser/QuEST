#include "aim_interaction.h"

AIM::AimInteraction::AimInteraction(const int interp_order,
                                    const Grid &grid,
                                    normalization::SpatialNorm normalization)
    : AimInteraction(nullptr,
                     nullptr,
                     nullptr,
                     interp_order,
                     1,
                     1,
                     grid,
                     Array2<Expansion>(),
                     normalization)
{
}

AIM::AimInteraction::AimInteraction(
    const std::shared_ptr<const DotVector> &dots,
    const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> &history,
    const std::shared_ptr<Propagation::RotatingFramePropagator> &propagator,
    const int interp_order,
    const double c0,
    const double dt,
    const Grid &grid,
    const Array2<Expansion> &expansion_table,
    normalization::SpatialNorm normalization)
    : HistoryInteraction(dots, history, propagator, interp_order, c0, dt),
      grid(grid),
      expansion_table(expansion_table),
      normalization(std::move(normalization)),
      max_transit_steps(grid.max_transit_steps(c0, dt)),
      circulant_dimensions(grid.circulant_shape(c0, dt, interp_order)),
      fourier_table(circulant_fourier_table()),
      source_table(circulant_dimensions),
      obs_table(circulant_dimensions),
      spatial_vector_transforms(spatial_fft_plans())
{
  std::fill(source_table.data(),
            source_table.data() + source_table.num_elements(),
            Eigen::Vector3cd::Zero());
  std::fill(obs_table.data(), obs_table.data() + obs_table.num_elements(),
            Eigen::Vector3cd::Zero());
}

const Interaction::ResultArray &AIM::AimInteraction::evaluate(const int step)
{
  // const auto wrapped_step = step % circulant_dimensions[0];
  // const auto nb = 8 * grid.num_boxes;

  // Eigen::Map<Eigen::ArrayXcd> observers(&obs_table[wrapped_step][0][0][0],
  // nb);
  // observers = 0;

  // for(int i = 1; i < circulant_dimensions[0]; ++i) {
  // if(step - i < 0) continue;
  // auto wrap = (step - i) % circulant_dimensions[0];

  // Eigen::Map<Eigen::ArrayXcd> prop(&fourier_table[i][0][0][0], nb);
  // Eigen::Map<Eigen::ArrayXcd> src(&source_table[wrap][0][0][0], nb);

  // observers += prop * src;
  //}

  // fftw_execute_dft(spatial_transforms.backward,
  // reinterpret_cast<fftw_complex *>(observers.data()),
  // reinterpret_cast<fftw_complex *>(observers.data()));

  // fill_results_table(step);
  return results;
}

void AIM::AimInteraction::fill_source_table(const int step)
{
  const int wrapped_step = step % circulant_dimensions[0];
  // std::cout << "(" << step << ", " << wrapped_step << ") ";
  auto p = &source_table[wrapped_step][0][0][0];
  std::fill(p, p + 8 * grid.num_boxes, Eigen::Vector3cd::Zero());

  for(auto dot_idx = 0u; dot_idx < expansion_table.shape()[0]; ++dot_idx) {
    for(auto expansion_idx = 0u; expansion_idx < expansion_table.shape()[1];
        ++expansion_idx) {
      const Expansion &e = expansion_table[dot_idx][expansion_idx];
      Eigen::Vector3i coord = grid.idx_to_coord(e.index);

      // This is the seam between what's stored in the History (density matrix
      // elements) and the electromagnetic source quantities. Ideally the AIM
      // code should not have knowledge of this to better encapsulate
      // "propagation," but this is good enough for now.
      source_table[wrapped_step][coord(0)][coord(1)][coord(2)] =
          e.weight * (*dots)[dot_idx].dipole() *
          history->array[dot_idx][step][0][RHO_01];
    }
  }

  fftw_execute_dft(
      spatial_vector_transforms.forward,
      reinterpret_cast<fftw_complex *>(&source_table[wrapped_step][0][0][0]),
      reinterpret_cast<fftw_complex *>(&source_table[wrapped_step][0][0][0]));
}

void AIM::AimInteraction::fill_results_table(const int step)
{
  results = 0;
  const int wrapped_step = step % circulant_dimensions[0];

  for(auto dot_idx = 0u; dot_idx < expansion_table.shape()[0]; ++dot_idx) {
    for(auto expansion_idx = 0u; expansion_idx < expansion_table.shape()[1];
        ++expansion_idx) {
      const Expansion &e = expansion_table[dot_idx][expansion_idx];
      Eigen::Vector3i coord = grid.idx_to_coord(e.index);

      results(dot_idx) +=
          e.weight *
          obs_table[wrapped_step][coord(0)][coord(1)][coord(2)].dot(
              (*dots)[dot_idx].dipole());
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
