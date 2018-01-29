#include <boost/test/unit_test.hpp>
#include <iostream>

#include "common.h"
#include "fourier.h"
#include "interactions/AIM/spacetime.h"

BOOST_AUTO_TEST_SUITE(SPACETIME)

struct SHAPE {
  std::array<int, 4> dims;
  SHAPE() : dims({{2, 6, 8, 10}}){};
};

BOOST_FIXTURE_TEST_CASE(CIRCULANT_MIRROR, SHAPE)
{
  spacetime::vector<int> stv(dims);
  std::fill(stv.data(), stv.data() + stv.num_elements(), 0);

  int i = 1;
  for(int t = 0; t < dims[0] / 2; ++t) {
    for(int x = 0; x < dims[1] / 2; ++x) {
      for(int y = 0; y < dims[2] / 2; ++y) {
        for(int z = 0; z < dims[3] / 2; ++z) {
          stv[t][x][y][z] = i++;
        }
      }
    }
  }

  spacetime::fill_circulant_mirror(stv);

  // Built using substitution rules in Mathematica -- does not contain the
  // flipped temporal entries!
  std::vector<int> check = {
      1,  2,  3,  4,  5,  0,  5,  4,  3,  2,  6,  7,  8,  9,  10, 0,  10, 9,
      8,  7,  11, 12, 13, 14, 15, 0,  15, 14, 13, 12, 16, 17, 18, 19, 20, 0,
      20, 19, 18, 17, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  16, 17, 18, 19,
      20, 0,  20, 19, 18, 17, 11, 12, 13, 14, 15, 0,  15, 14, 13, 12, 6,  7,
      8,  9,  10, 0,  10, 9,  8,  7,  21, 22, 23, 24, 25, 0,  25, 24, 23, 22,
      26, 27, 28, 29, 30, 0,  30, 29, 28, 27, 31, 32, 33, 34, 35, 0,  35, 34,
      33, 32, 36, 37, 38, 39, 40, 0,  40, 39, 38, 37, 0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  36, 37, 38, 39, 40, 0,  40, 39, 38, 37, 31, 32, 33, 34,
      35, 0,  35, 34, 33, 32, 26, 27, 28, 29, 30, 0,  30, 29, 28, 27, 41, 42,
      43, 44, 45, 0,  45, 44, 43, 42, 46, 47, 48, 49, 50, 0,  50, 49, 48, 47,
      51, 52, 53, 54, 55, 0,  55, 54, 53, 52, 56, 57, 58, 59, 60, 0,  60, 59,
      58, 57, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  56, 57, 58, 59, 60, 0,
      60, 59, 58, 57, 51, 52, 53, 54, 55, 0,  55, 54, 53, 52, 46, 47, 48, 49,
      50, 0,  50, 49, 48, 47, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  41, 42, 43, 44,
      45, 0,  45, 44, 43, 42, 46, 47, 48, 49, 50, 0,  50, 49, 48, 47, 51, 52,
      53, 54, 55, 0,  55, 54, 53, 52, 56, 57, 58, 59, 60, 0,  60, 59, 58, 57,
      0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  56, 57, 58, 59, 60, 0,  60, 59,
      58, 57, 51, 52, 53, 54, 55, 0,  55, 54, 53, 52, 46, 47, 48, 49, 50, 0,
      50, 49, 48, 47, 21, 22, 23, 24, 25, 0,  25, 24, 23, 22, 26, 27, 28, 29,
      30, 0,  30, 29, 28, 27, 31, 32, 33, 34, 35, 0,  35, 34, 33, 32, 36, 37,
      38, 39, 40, 0,  40, 39, 38, 37, 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
      36, 37, 38, 39, 40, 0,  40, 39, 38, 37, 31, 32, 33, 34, 35, 0,  35, 34,
      33, 32, 26, 27, 28, 29, 30, 0,  30, 29, 28, 27};

  BOOST_CHECK_EQUAL(check.size(), stv.num_elements() / 2);
  for(size_t i = 0; i < check.size(); ++i) {
    BOOST_CHECK_EQUAL(check.at(i), stv.data()[i]);
  }
}

BOOST_FIXTURE_TEST_CASE(FFT_MATRIX_MULTIPLY, SHAPE)
{
  const int num_spatial_elements = dims[1] * dims[2] * dims[3];
  spacetime::vector<cmplx> mat(dims), vec(dims);

  TransformPair transforms = {
      fftw_plan_dft_3d(dims[1], dims[2], dims[3],
                       reinterpret_cast<fftw_complex *>(mat.data()),
                       reinterpret_cast<fftw_complex *>(mat.data()),
                       FFTW_FORWARD, FFTW_MEASURE),
      fftw_plan_dft_3d(dims[1], dims[2], dims[3],
                       reinterpret_cast<fftw_complex *>(mat.data()),
                       reinterpret_cast<fftw_complex *>(mat.data()),
                       FFTW_BACKWARD, FFTW_MEASURE)};

  std::fill(mat.data(), mat.data() + mat.num_elements(), cmplx(0, 0));
  std::fill(vec.data(), vec.data() + vec.num_elements(), cmplx(0, 0));

  int i = 1;
  for(int t = 0; t < dims[0] / 2; ++t) {
    for(int x = 0; x < dims[1] / 2; ++x) {
      for(int y = 0; y < dims[2] / 2; ++y) {
        for(int z = 0; z < dims[3] / 2; ++z) {
          mat[0][x][y][z] = i;
          vec[0][x][y][z] = i;
          ++i;
        }
      }
    }
  }

  spacetime::fill_circulant_mirror(mat);

  fftw_execute_dft(transforms.forward,
                   reinterpret_cast<fftw_complex *>(mat.data()),
                   reinterpret_cast<fftw_complex *>(mat.data()));
  fftw_execute_dft(transforms.forward,
                   reinterpret_cast<fftw_complex *>(vec.data()),
                   reinterpret_cast<fftw_complex *>(vec.data()));

  spacetime::vector<cmplx> result(dims);
  Eigen::Map<Eigen::ArrayXcd> fft_mat(mat.data(), num_spatial_elements);
  Eigen::Map<Eigen::ArrayXcd> fft_vec(vec.data(), num_spatial_elements);
  Eigen::Map<Eigen::ArrayXcd> fft_result(result.data(), num_spatial_elements);
  fft_mat /= num_spatial_elements;
  fft_result = fft_mat * fft_vec;

  fftw_execute_dft(transforms.backward,
                   reinterpret_cast<fftw_complex *>(result.data()),
                   reinterpret_cast<fftw_complex *>(result.data()));

  std::vector<int> check = {
      73810, 72664, 72226, 72520, 73570, 68110, 66964, 66526, 66820, 67870,
      66610, 65464, 65026, 65320, 66370, 70060, 68914, 68476, 68770, 69820,
      45610, 44464, 44026, 44320, 45370, 39910, 38764, 38326, 38620, 39670,
      38410, 37264, 36826, 37120, 38170, 41860, 40714, 40276, 40570, 41620,
      41810, 40664, 40226, 40520, 41570, 36110, 34964, 34526, 34820, 35870,
      34610, 33464, 33026, 33320, 34370, 38060, 36914, 36476, 36770, 37820};

  const double TOLER = 1e-13;
  int j = 0;
  for(int t = 0; t < dims[0] / 2; ++t) {
    for(int x = 0; x < dims[1] / 2; ++x) {
      for(int y = 0; y < dims[2] / 2; ++y) {
        for(int z = 0; z < dims[3] / 2; ++z) {
          BOOST_CHECK_CLOSE(result[t][x][y][z].real(), check.at(j++), TOLER);
          BOOST_CHECK_CLOSE(result[t][x][y][z].imag(), 0, TOLER);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()  // common
