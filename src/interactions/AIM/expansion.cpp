#include "expansion.h"

Eigen::MatrixXd AIM::Expansions::LeastSquaresExpansionSolver::w_matrix(
    const Eigen::Vector3d &pos) const
{
  Eigen::MatrixXd w_mat = Eigen::MatrixXd::Zero(num_pts, num_pts);

  auto expansion_indices = grid.expansion_indices(pos);

  for(int col = 0; col < num_pts; ++col) {
    Eigen::Vector3d dr =
        grid.spatial_coord_of_box(expansion_indices.at(col)) - pos;
    int row = 0;
    for(int nx = 0; nx <= grid.order(); ++nx) {
      double x_term = std::pow(dr(0), nx);
      for(int ny = 0; ny <= grid.order(); ++ny) {
        double y_term = std::pow(dr(1), ny);
        for(int nz = 0; nz <= grid.order(); ++nz) {
          double z_term = std::pow(dr(2), nz);
          w_mat(row++, col) = x_term * y_term * z_term;
        }
      }
    }
  }

  return w_mat;
}

Eigen::VectorXd AIM::Expansions::LeastSquaresExpansionSolver::solve(
    const Eigen::Vector3d &pos) const
{
  Eigen::FullPivLU<Eigen::MatrixXd> lu(w_matrix(pos));
  return lu.solve(qvec);
}

AIM::Expansions::ExpansionTable
AIM::Expansions::LeastSquaresExpansionSolver::table(
    const std::vector<QuantumDot> &dots) const
{
  AIM::Expansions::ExpansionTable table(boost::extents[dots.size()][num_pts]);

  for(int dot_idx = 0; dot_idx < static_cast<int>(dots.size()); ++dot_idx) {
    const auto weights = solve(dots.at(dot_idx).position());
    const auto indices = grid.expansion_indices(dots.at(dot_idx).position());

    for(auto w = 0; w < num_pts; ++w)
      table[dot_idx][w] = {indices[w], weights(w)};
  }

  return table;
}

boost::multi_array<double, 4>
AIM::Expansions::LeastSquaresExpansionSolver::chebyshev_table(
    const std::vector<double> &pts) const
{
  boost::multi_array<double, 4> coefs(
      boost::extents[pts.size()][pts.size()][pts.size()][num_pts]);

  const int n = pts.size();
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      for(int k = 0; k < n; ++k) {
        Eigen::Vector3d r =
            Eigen::Array3d(pts[i], pts[j], pts[k]) * grid.spacing();
        const auto weights = solve(r);
        for(int w = 0; w < num_pts; ++w) coefs[i][j][k][w] = weights[w];
      }
    }
  }

  return coefs;
}
