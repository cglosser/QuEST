#ifndef BLOCH_SYSTEM_H
#define BLOCH_SYSTEM_H

#include <utility>
#include <vector>

#include "integrator.h"
#include "interaction_table.h"
#include "quantum_dot.h"

class BlochSystem {
 public:
  BlochSystem(const PredictorCorrector, std::vector<QuantumDot>, const size_t);
  void step();

  // private:
  typedef boost::multi_array<matrix_elements, 3> HistoryArray;

  PredictorCorrector integrator;
  std::vector<QuantumDot> dots;
  HistoryArray history;
  const int num_dots;
  const int num_steps;
};
#endif
