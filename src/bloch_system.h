#ifndef BLOCH_SYSTEM_H
#define BLOCH_SYSTEM_H

#include <Eigen/Dense>
#include <utility>
#include <vector>

#include "integrator.h"
#include "interaction_table.h"
#include "quantum_dot.h"

class BlochSystem {
 public:
  BlochSystem(const PredictorCorrector::Weights, std::vector<QuantumDot>,
              const int, const size_t);
  void step();

  // private:
  typedef boost::multi_array<matrix_elements, 3,
                             Eigen::aligned_allocator<matrix_elements>>
      HistoryArray;

  const int num_steps;
  const int num_dots;

  PredictorCorrector::Weights integrator;
  std::vector<QuantumDot> dots;
  HistoryArray history;
  InteractionTable interactions;

  std::vector<double> rabi_freqs;

  int now;
  const double dt;

  void convolve_currents(const int time_idx);
};
#endif
