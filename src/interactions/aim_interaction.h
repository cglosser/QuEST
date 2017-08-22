#ifndef AIM_INTERACTION_H
#define AIM_INTERACTION_H

#include <fftw3.h>
#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>

#include <iomanip>

#include "../common.h"
#include "history_interaction.h"

namespace AIM {
  class Grid;
  class AimInteraction;
}

class AIM::Grid {
 public:
  typedef Eigen::Array<int, 3, 2> BoundsArray;

  Grid();
  Grid(const Eigen::Array3d &, const std::shared_ptr<DotVector> &, const int);
  Grid(const Eigen::Array3d &, const Eigen::Array3i &);

  std::vector<DotRange> box_contents_map(
      const std::shared_ptr<DotVector> &) const;

  // Geometry routines (grid <---> space)
  Eigen::Vector3i grid_coordinate(const Eigen::Vector3d &) const;
  size_t coord_to_idx(const Eigen::Vector3i &) const;
  Eigen::Vector3i idx_to_coord(size_t) const;
  Eigen::Vector3d spatial_coord_of_box(const size_t) const;
  std::vector<size_t> expansion_box_indices(const Eigen::Vector3d &,
                                            const int) const;

  Eigen::Array3i dimensions;
  size_t num_boxes;
  double max_diagonal;

  int max_transit_steps(double c, double dt) const
  {
    return static_cast<int>(ceil(max_diagonal / (c * dt)));
  };

  BoundsArray calculate_bounds() const;

 private:
  Eigen::Array3d spacing;
  std::shared_ptr<DotVector> dots;
  int expansion_order;
  BoundsArray bounds;

  void sort_points_on_boxidx() const;
};

class AIM::AimInteraction final : public HistoryInteraction {
 public:
  struct Expansion {
    size_t index;
    double weight;
  };

  AimInteraction(
      const std::shared_ptr<const DotVector> &,
      const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> &,
      const std::shared_ptr<Propagation::RotatingFramePropagator> &,
      const int,
      const double,
      const double,
      const Grid &,
      const int);

  const ResultArray &evaluate(const int) final;
  void fill_gmatrix_table(SpacetimeVector<double> &) const;
  Eigen::VectorXd q_vector(const Eigen::Vector3d &) const;
  Eigen::MatrixXd w_matrix(const Eigen::Vector3d &) const;
  Eigen::VectorXd solve_expansion_system(const Eigen::Vector3d &) const;

  Array<Expansion> expansions() const;

 //private:
  Grid grid;
  int box_order, max_transit_steps;

  SpacetimeVector<cmplx> fourier_table;
  Array<Expansion> expansion_table;
  boost::multi_array<cmplx, 2> source_table;
  fftw_plan vector_forward_plan, vector_backward_plan;

  void fill_fourier_table();
  std::pair<fftw_plan, fftw_plan> vector_fft_plans();
  void fill_source_table(const int);
};

#endif
