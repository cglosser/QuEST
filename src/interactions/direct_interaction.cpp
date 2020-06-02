#include "direct_interaction.h"

DirectInteraction::DirectInteraction(
    std::shared_ptr<const DotVector> dots,
    std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history,
    Propagation::Kernel<cmplx> &kernel,
    const int interp_order,
    const double c0,
    const double dt)
    : HistoryInteraction(
          std::move(dots), std::move(history), interp_order, c0, dt),
      num_interactions((this->dots)->size() * ((this->dots)->size() - 1) / 2),
      floor_delays(num_interactions),
      coefficients(boost::extents[num_interactions][interp_order + 1])
{
  build_coefficient_table(kernel);
}

void DirectInteraction::build_coefficient_table(
    Propagation::Kernel<cmplx> &kernel)
{
  Interpolation::UniformLagrangeSet lagrange(interp_order);

  for(int pair_idx = 0; pair_idx < num_interactions; ++pair_idx) {
    int src, obs;
    std::tie(src, obs) = idx2coord(pair_idx);

    Eigen::Vector3d dr(separation((*dots)[src], (*dots)[obs]));

    std::pair<int, double> delay(split_double(dr.norm() / (c0 * dt)));

    floor_delays[pair_idx] = delay.first;
    lagrange.evaluate_derivative_table_at_x(delay.second, dt);

    std::vector<Eigen::Matrix3cd> interp_dyads(
        kernel.coefficients(dr, lagrange));

    for(int i = 0; i <= interp_order; ++i) {
      coefficients[pair_idx][i] =
          (*dots)[obs].dipole().dot(interp_dyads[i] * (*dots)[src].dipole());
    }
  }
}

const InteractionBase::ResultArray &
DirectInteraction::first_evaluation_of_timestep(const int time_idx)
{
  constexpr int RHO_01 = 1;
  results.setZero();
  past_terms_of_results.setZero();

  // iterate through all dot pairs
  for(int pair_idx = 0; pair_idx < num_interactions; ++pair_idx) {
    int src, obs;
    std::tie(src, obs) = idx2coord(pair_idx);

    for(int i = 1; i <= interp_order; ++i) {
      const int s =
          std::max(time_idx - floor_delays[pair_idx] - i,
                   static_cast<int>(history->array_.index_bases()[1]));

      past_terms_of_results[src] +=
          (history->array_[obs][s][0])[RHO_01] * coefficients[pair_idx][i];
      past_terms_of_results[obs] +=
          (history->array_[src][s][0])[RHO_01] * coefficients[pair_idx][i];
    }
    const int s = std::max(time_idx - floor_delays[pair_idx],
                           static_cast<int>(history->array_.index_bases()[1]));

    results[src] +=
        (history->array_[obs][s][0])[RHO_01] * coefficients[pair_idx][0];
    results[obs] +=
        (history->array_[src][s][0])[RHO_01] * coefficients[pair_idx][0];
  }
  results += past_terms_of_results;
  return results;
}

const InteractionBase::ResultArray &DirectInteraction::evaluate(
    const int time_idx)
{
  constexpr int RHO_01 = 1;
  results.setZero();

  // iterate through all dot pairs
  for(int pair_idx = 0; pair_idx < num_interactions; ++pair_idx) {
    int src, obs;
    std::tie(src, obs) = idx2coord(pair_idx);

    const int s = std::max(time_idx - floor_delays[pair_idx],
                           static_cast<int>(history->array_.index_bases()[1]));

    results[src] +=
        (history->array_[obs][s][0])[RHO_01] * coefficients[pair_idx][0];
    results[obs] +=
        (history->array_[src][s][0])[RHO_01] * coefficients[pair_idx][0];
  }
  results += past_terms_of_results;
  return results;
}

// const InteractionBase::ResultArray &DirectInteraction::evaluate(
//     const int time_idx, const bool first_call)
// {
//   constexpr int RHO_01 = 1;
//   temp_res.setZero();
//   if(first_call) results.setZero();
//
//   // iterate through all dot pairs
//   for(int pair_idx = 0; pair_idx < num_interactions; ++pair_idx) {
//     int src, obs;
//     std::tie(src, obs) = idx2coord(pair_idx);
//
//     if(first_call) {  // sum over history if this is the first call
//       for(int i = 1; i <= interp_order; ++i) {
//         const int s =
//             std::max(time_idx - floor_delays[pair_idx] - i,
//                      static_cast<int>(history->array_.index_bases()[1]));
//
//         results[src] +=
//             (history->array_[obs][s][0])[RHO_01] * coefficients[pair_idx][i];
//         results[obs] +=
//             (history->array_[src][s][0])[RHO_01] * coefficients[pair_idx][i];
//       }
//     }
//     const int s = std::max(time_idx - floor_delays[pair_idx],
//                            static_cast<int>(history->array_.index_bases()[1]));
//
//     temp_res[src] +=
//         (history->array_[obs][s][0])[RHO_01] * coefficients[pair_idx][0];
//     temp_res[obs] +=
//         (history->array_[src][s][0])[RHO_01] * coefficients[pair_idx][0];
//   }
//   temp_res += results;
//   return temp_res;
// }

int DirectInteraction::coord2idx(int row, int col)
{
  assert(row != col);
  if(col > row) std::swap(row, col);

  return row * (row - 1) / 2 + col;
}

std::pair<int, int> DirectInteraction::idx2coord(const int idx)
{
  const int row = std::floor((std::sqrt(1 + 8 * idx) + 1) / 2.0);
  const int col = idx - row * (row - 1) / 2;

  return std::pair<int, int>(row, col);
}
