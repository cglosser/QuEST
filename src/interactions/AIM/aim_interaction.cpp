#include "aim_interaction.h"

std::vector<DotRange> AIM::Grid::box_contents_map(
    const std::shared_ptr<DotVector> &dots) const
{
  std::vector<DotRange> boxes(num_boxes);
  for(size_t box_idx = 0; box_idx < boxes.size(); ++box_idx) {
    auto IsInBox = [=](const QuantumDot &qd) {
      return coord_to_idx(grid_coordinate(qd.position())) == box_idx;
    };

    auto begin = std::find_if(dots->begin(), dots->end(), IsInBox);
    auto end = std::find_if_not(begin, dots->end(), IsInBox);
    boxes.at(box_idx) = std::make_pair(begin, end);
  }

  return boxes;
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

AIM::AimInteraction::AimInteraction(
    const std::shared_ptr<const DotVector> &dots,
    const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> &history,
    const std::shared_ptr<Propagation::RotatingFramePropagator> &propagator,
    const int interp_order,
    const double c0,
    const double dt,
    const Grid &grid,
    const int box_order)
    : HistoryInteraction(dots, history, propagator, interp_order, c0, dt),
      grid(grid),
      box_order(box_order),
      max_transit_steps(grid.max_transit_steps(c0, dt)),
      circulant_dimensions(grid.circulant_shape(c0, dt)),
      expansion_table(expansions()),
      fourier_table(circulant_fourier_table()),
      source_table(circulant_dimensions),
      obs_table(circulant_dimensions),
      spatial_transforms(spatial_fft_plans())
{
  // std::fill(source_table.data(),
  // source_table.data() + source_table.num_elements(), cmplx(0, 0));
  // std::fill(obs_table.data(), obs_table.data() + obs_table.num_elements(),
  // cmplx(0, 0));
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
  for(size_t dot_idx = 0; dot_idx < dots->size(); ++dot_idx) {
    for(int idx = 0; idx < std::pow(box_order + 1, 3); ++idx) {
      const Expansion &e = expansion_table[dot_idx][idx];
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
  fftw_plan circulant_plan = fftw_plan_many_dft(
      3, &circulant_dimensions[1], circulant_dimensions[0],
      reinterpret_cast<fftw_complex *>(g_mat.data()), NULL, 1, num_gridpts,
      reinterpret_cast<fftw_complex *>(g_mat.data()), NULL, 1, num_gridpts,
      FFTW_FORWARD, FFTW_MEASURE);

  std::fill(g_mat.data(), g_mat.data() + g_mat.num_elements(), cmplx(0, 0));

  fill_gmatrix_table(g_mat);

  // Transform the circulant vectors into their equivalently-diagonal
  // representation. Buckle up.

  fftw_execute(circulant_plan);

  // This accounts for FFTW's *un*normalized transform -- it takes the least
  // amount of computational effort to put all of the normalizations here.

  Eigen::Map<Eigen::ArrayXcd> gs(g_mat.data(), g_mat.num_elements());
  gs /= num_gridpts;

  fftw_destroy_plan(circulant_plan);

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
            gmatrix_table[t][x][y][z] = interp.evaluations[0][polynomial_idx];
          }
        }
      }
    }
  }

  fill_circulant_mirror(gmatrix_table);
}

AIM::AimInteraction::TransformPair AIM::AimInteraction::spatial_fft_plans()
{
  // Set up FFTW plans to transform projected source distributions. Due to the
  // requirements of the circulant extension, these plans perform transforms of
  // length 2 n_{x,y,z} to accommodate the requisite zero padding. While they're
  // constructed to work on the head of `source_table` (that is, what would be
  // the I_0 source), the advanced FFTW interface allows them to stride forward
  // to equivalently transform the source currents at every timestep.

  auto make_plan = [=](const int sign) {
    return fftw_plan_dft_3d(
        circulant_dimensions[1], circulant_dimensions[2],
        circulant_dimensions[3],
        reinterpret_cast<fftw_complex *>(source_table.data()),
        reinterpret_cast<fftw_complex *>(source_table.data()), sign,
        FFTW_MEASURE);
  };

  fftw_plan fwd, bkwd;

  fwd = make_plan(FFTW_FORWARD);
  bkwd = make_plan(FFTW_BACKWARD);

  return {fwd, bkwd};
}

Array<AIM::AimInteraction::Expansion> AIM::AimInteraction::expansions() const
{
  const size_t num_pts = std::pow(box_order + 1, 3);
  Array<AIM::AimInteraction::Expansion> table(
      boost::extents[dots->size()][num_pts]);

  for(size_t dot_idx = 0; dot_idx < dots->size(); ++dot_idx) {
    const auto indices =
        grid.expansion_box_indices(dots->at(dot_idx).position(), box_order);
    const auto weights = solve_expansion_system(dots->at(dot_idx).position());

    for(size_t w = 0; w < num_pts; ++w) {
      table[dot_idx][w] = {indices[w], weights[w]};
    }
  }

  return table;
}

Eigen::VectorXd AIM::AimInteraction::q_vector(const Eigen::Vector3d &pos) const
{
  const int len = std::pow(box_order + 1, 3);
  Eigen::VectorXd q_vec(len);

  int i = 0;
  for(int nx = 0; nx <= box_order; ++nx) {
    for(int ny = 0; ny <= box_order; ++ny) {
      for(int nz = 0; nz <= box_order; ++nz) {
        q_vec(i++) =
            std::pow(pos(0), nx) * std::pow(pos(1), ny) * std::pow(pos(2), nz);
      }
    }
  }

  return q_vec;
}

Eigen::MatrixXd AIM::AimInteraction::w_matrix(const Eigen::Vector3d &pos) const
{
  const int len = std::pow(box_order + 1, 3);
  Eigen::MatrixXd w_mat(len, len);

  auto expansion_indices = grid.expansion_box_indices(pos, box_order);

  for(int i = 0; i < len; ++i) {
    w_mat.col(i) = q_vector(grid.spatial_coord_of_box(expansion_indices.at(i)));
  }

  return w_mat;
}

Eigen::VectorXd AIM::AimInteraction::solve_expansion_system(
    const Eigen::Vector3d &pos) const
{
  Eigen::FullPivLU<Eigen::MatrixXd> lu(w_matrix(pos));
  return lu.solve(q_vector(pos));
}
