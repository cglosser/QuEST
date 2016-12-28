#ifndef PULSE_INTERACTION_H
#define PULSE_INTERACTION_H

#include "../configuration.h"
#include "../pulse.h"
#include "interaction.h"

class PulseInteraction : public Interaction {
 public:
  PulseInteraction(const std::shared_ptr<const DotVector> &,
                   const std::shared_ptr<const Pulse>);
  void evaluate(const int);

 private:
  std::shared_ptr<const Pulse> pulse;
};

#endif
