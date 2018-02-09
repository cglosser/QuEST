#ifndef SPACETIME_H
#define SPACETIME_H

namespace spacetime {
  template <class T>
  using vector = boost::multi_array<T, 4>;

  template <class T>
  using vector3d = boost::multi_array<T, 5>;
  // Last dim holds X, Y, Z components of a vector quantity

  template <class T>
  vector3d<T> make_vector3d(const std::array<int, 4> &dims)
  {
    // Convenience function to build vector3d from "circulant shape" array
    return vector3d<T>(boost::extents[dims[0]][dims[1]][dims[2]][dims[3]][3]);
  }

  template <class T>
  void fill_circulant_mirror(spacetime::vector<T> &);

  template <class T>
  void clear_time_slice(spacetime::vector3d<T> &, const int);
}

template <class T>
void spacetime::fill_circulant_mirror(spacetime::vector<T> &stv)
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

template <class T>
void spacetime::clear_time_slice(spacetime::vector3d<T> &table, const int s)
{
  const auto ptr = &table[s][0][0][0][0];
  const auto n_elem = std::accumulate(table.shape() + 1, table.shape() + 5, 1,
                                      std::multiplies<int>());
  std::fill(ptr, ptr + n_elem, cmplx(0, 0));
}

#endif
