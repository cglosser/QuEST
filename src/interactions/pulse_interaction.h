#ifndef PULSE_INTERACTION_H
#define PULSE_INTERACTION_H

#include "interaction.h"

class PulseInteraction : Interaction {
 public:
  PulseInteraction(const std::shared_ptr<const DotVector> &,
                   const std::shared_ptr<const Pulse>);
  void evaluate(const int);

 private:
  std::shared_ptr<const Pulse> pulse;
};

#endif
