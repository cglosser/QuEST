#ifndef BLOCH_SYSTEM_H
#define BLOCH_SYSTEM_H

#include <utility>
#include <vector>

#include "integrator.h"
#include "interaction_table.h"
#include "quantum_dot.h"

class BlochSystem {
 public:
  BlochSystem(std::vector<QuantumDot>);
  void step();

 //private:
  std::vector<QuantumDot> dots;
};
#endif
