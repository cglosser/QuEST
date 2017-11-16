#ifndef AIM_INTERACTION_H
#define AIM_INTERACTION_H

#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>

#include "common.h"
#include "expansion.h"
#include "fourier.h"
#include "grid.h"
#include "interactions/history_interaction.h"

namespace AIM {
  class AimInteraction;

  namespace spacetime {
    template <class T>
    using vector = boost::multi_array<T, 4>;

    template <class T>
    using vector3d = boost::multi_array<T, 5>;
    // Last dim holds X, Y, Z components of a vector quantity

    template <class T>
    void fill_circulant_mirror(spacetime::vector<T> &);
  }

  namespace normalization {
    using SpatialNorm = std::function<double(const Eigen::Vector3d &)>;
    const SpatialNorm unit = [](__attribute__((unused))
                                const Eigen::Vector3d &v) { return 1; };
    const SpatialNorm distance = [](const Eigen::Vector3d &v) {
      return v.norm();
    };
    const SpatialNorm poisson = [](const Eigen::Vector3d &v) {
      return 4 * M_PI * v.norm();
    };
  }
}

class AIM::AimInteraction final : public HistoryInteraction {
 public:
  AimInteraction(const int, const Grid &, normalization::SpatialNorm);
  AimInteraction(
      const std::shared_ptr<const DotVector> &,
      const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> &,
      const std::shared_ptr<Propagation::RotatingFramePropagator> &,
      const int,
      const double,
      const double,
      const Grid &,
      const Expansions::ExpansionTable &,
      normalization::SpatialNorm);

  const ResultArray &evaluate(const int) final;

  // private:
  Grid grid;
  Expansions::ExpansionTable expansion_table;
  normalization::SpatialNorm normalization;
  int max_transit_steps;
  std::array<int, 4> circulant_dimensions;

  // This corresponds to delta(t - R/c)/R and thus holds *scalar* quantities
  spacetime::vector<cmplx> fourier_table;

  // These correspond to J and E and thus hold *vector* quantities
  spacetime::vector<Eigen::Vector3cd> source_table, obs_table;

  TransformPair spatial_vector_transforms;

  void fill_source_table(const int);
  void fill_results_table(const int);

  spacetime::vector<cmplx> circulant_fourier_table();
  void fill_gmatrix_table(spacetime::vector<cmplx> &) const;
  TransformPair spatial_fft_plans();
};

template <class T>
void AIM::spacetime::fill_circulant_mirror(spacetime::vector<T> &stv)
{
  // This routine does *NOT* mirror the temporal "axis" (on purpose)
  // ...it's also really, really, clever ;)
  const int nx = stv.shape()[1] / 2;
  const int ny = stv.shape()[2] / 2;
  const int nz = stv.shape()[3] / 2;

  const int x_stride = (2 * ny) * (2 * nz);
  const int y_stride = (2 * nz);
  const int z_stride = 1;

  using MapArray = Eigen::Map<Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>>;

  const int start = static_cast<int>(stv.index_bases()[0]);
  const int stop = static_cast<int>(start + stv.shape()[0]);

  for(int t = start; t < stop; ++t) {
    for(int x = 0; x < nx; ++x) {
      for(int y = 0; y < ny; ++y) {
        MapArray zflip(&stv[t][x][y][1], z_stride, nz - 1);
        MapArray zhole(&stv[t][x][y][nz + 1], z_stride, nz - 1);
        zhole = zflip.rowwise().reverse();
      }
      MapArray yflip(&stv[t][x][1][0], y_stride, ny - 1);
      MapArray yhole(&stv[t][x][ny + 1][0], y_stride, ny - 1);
      yhole = yflip.rowwise().reverse();
    }
    MapArray xflip(&stv[t][1][0][0], x_stride, nx - 1);
    MapArray xhole(&stv[t][nx + 1][0][0], x_stride, nx - 1);
    xhole = xflip.rowwise().reverse();
  }
}

#endif
