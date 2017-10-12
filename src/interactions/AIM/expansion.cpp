#include "expansion.h"

Array<AIM::Expansion> AIM::LeastSquaresExpansionSolver::get(
    const int box_order, const Grid &grid, const std::vector<QuantumDot> &dots)
{
  return LeastSquaresExpansionSolver(box_order, grid).table(dots);
}

Array<AIM::Expansion> AIM::LeastSquaresExpansionSolver::table(
    const std::vector<QuantumDot> &dots) const
{
  Array<AIM::Expansion> table(boost::extents[dots.size()][num_pts]);

  for(auto dot_idx = 0u; dot_idx < dots.size(); ++dot_idx) {
    const auto indices =
        grid.expansion_box_indices(dots.at(dot_idx).position(), box_order);
    const auto weights = solve_expansion_system(dots.at(dot_idx).position());

    for(auto w = 0; w < num_pts; ++w) {
      table[dot_idx][w] = {indices[w], weights[w]};
    }
  }

  return table;
}

Eigen::VectorXd AIM::LeastSquaresExpansionSolver::q_vector(
    const Eigen::Vector3d &pos) const
{
  Eigen::VectorXd q_vec(num_pts);

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

Eigen::MatrixXd AIM::LeastSquaresExpansionSolver::w_matrix(
    const Eigen::Vector3d &pos) const
{
  Eigen::MatrixXd w_mat(num_pts, num_pts);

  auto expansion_indices = grid.expansion_box_indices(pos, box_order);

  for(int i = 0; i < num_pts; ++i) {
    w_mat.col(i) = q_vector(grid.spatial_coord_of_box(expansion_indices.at(i)));
  }

  return w_mat;
}

Eigen::VectorXd AIM::LeastSquaresExpansionSolver::solve_expansion_system(
    const Eigen::Vector3d &pos) const
{
  Eigen::FullPivLU<Eigen::MatrixXd> lu(w_matrix(pos));
  return lu.solve(q_vector(pos));
}
