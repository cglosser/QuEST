#include "expansion.h"

AIM::Expansions::ExpansionTable
AIM::Expansions::LeastSquaresExpansionSolver::get_expansions(
    const int box_order, const Grid &grid, const std::vector<QuantumDot> &dots)
{
  return LeastSquaresExpansionSolver(box_order, grid).table(dots);
}

AIM::Expansions::ExpansionTable
AIM::Expansions::LeastSquaresExpansionSolver::table(
    const std::vector<QuantumDot> &dots) const
{
  using namespace enums;
  AIM::Expansions::ExpansionTable table(
      boost::extents[dots.size()][num_pts]);

  for(auto dot_idx = 0u; dot_idx < dots.size(); ++dot_idx) {
    const auto &pos = dots.at(dot_idx).position();
    const auto indices = grid.expansion_box_indices(pos, box_order);

    Eigen::Array4Xd weights(4, num_pts);
    weights.row(D_0) = solve_expansion_system(pos, {{0, 0, 0}});
    weights.row(D_X) = solve_expansion_system(pos, {{1, 0, 0}});
    weights.row(D_Y) = solve_expansion_system(pos, {{0, 1, 0}});
    weights.row(D_Z) = solve_expansion_system(pos, {{0, 0, 1}});

    for(auto w = 0; w < num_pts; ++w) {
      table[dot_idx][w].index = indices[w];
      Eigen::Map<Eigen::Array4d>(table[dot_idx][w].weights.data()) =
          weights.col(w);
    }
  }

  return table;
}

Eigen::VectorXd AIM::Expansions::LeastSquaresExpansionSolver::q_vector(
    const Eigen::Vector3d &pos, const std::array<int, 3> &derivatives) const
{
  Eigen::VectorXd q_vec(num_pts);

  int i = 0;
  for(int nx = 0; nx <= box_order; ++nx) {
    double x_term = falling_factorial(nx, derivatives[0]) *
                    std::pow(pos(0), nx - derivatives[0]);
    for(int ny = 0; ny <= box_order; ++ny) {
      double y_term = falling_factorial(ny, derivatives[1]) *
                      std::pow(pos(1), ny - derivatives[1]);
      for(int nz = 0; nz <= box_order; ++nz) {
        double z_term = falling_factorial(nz, derivatives[2]) *
                        std::pow(pos(2), nz - derivatives[2]);
        q_vec(i++) = x_term * y_term * z_term;
      }
    }
  }

  return q_vec;
}

Eigen::MatrixXd AIM::Expansions::LeastSquaresExpansionSolver::w_matrix(
    const Eigen::Vector3d &pos) const
{
  Eigen::MatrixXd w_mat(num_pts, num_pts);

  auto expansion_indices = grid.expansion_box_indices(pos, box_order);

  for(int i = 0; i < num_pts; ++i) {
    w_mat.col(i) = q_vector(grid.spatial_coord_of_box(expansion_indices.at(i)));
  }

  return w_mat;
}

Eigen::VectorXd
AIM::Expansions::LeastSquaresExpansionSolver::solve_expansion_system(
    const Eigen::Vector3d &pos, const std::array<int, 3> &derivatives) const
{
  Eigen::FullPivLU<Eigen::MatrixXd> lu(w_matrix(pos));
  return lu.solve(q_vector(pos, derivatives));
}
