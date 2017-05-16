#include "../src/quantum_dot.h"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(quantum_dot)

BOOST_AUTO_TEST_SUITE(liouville_rhs)

const Eigen::Vector3d pos(0, 0, 0), dipole(0, 0, 1);
const double omega = 2278.9013;
const std::pair<double, double> ts(10, 20);
const QuantumDot qd(pos, omega, ts, dipole);

BOOST_AUTO_TEST_CASE(resonant_rhs)
{
  const Eigen::Vector3cd efield(0, 0, 1);
  const cmplx rabi = efield.dot(dipole);

  Eigen::Matrix2cd hamiltonian_matrix;
  hamiltonian_matrix << 0, rabi, std::conj(rabi), 0;

  Eigen::Matrix2cd density_matrix;
  density_matrix << 1, cmplx(0.5, 0.5), cmplx(0.5, -0.5), 0;

  Eigen::Matrix2cd dissipator;
  dissipator << (density_matrix(0, 0) - 1.0) / ts.first,
      density_matrix(0, 1) / ts.second, density_matrix(1, 0) / ts.second,
      density_matrix(1, 1) / ts.first;

  // Evaluate the commutator & dissipator
  Eigen::Matrix2cd matrix_rhs(
      -iu * (hamiltonian_matrix * density_matrix -
             density_matrix * hamiltonian_matrix) - dissipator);

  matrix_elements rhs(qd.liouville_rhs(
      matrix_elements(density_matrix(0, 0), density_matrix(0, 1)), rabi,
      omega));

  BOOST_CHECK_CLOSE(matrix_rhs(0, 0).real(), rhs(0).real(), 1e-9);
  BOOST_CHECK_CLOSE(matrix_rhs(0, 0).imag(), rhs(0).imag(), 1e-9);
                                                                
  BOOST_CHECK_CLOSE(matrix_rhs(0, 1).real(), rhs(1).real(), 1e-9);
  BOOST_CHECK_CLOSE(matrix_rhs(0, 1).imag(), rhs(1).imag(), 1e-9);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
