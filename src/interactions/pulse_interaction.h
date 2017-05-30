#ifndef PULSE_INTERACTION_H
#define PULSE_INTERACTION_H

#include "../configuration.h"
#include "../pulse.h"
#include "interaction.h"

class PulseInteraction : public Interaction {
 public:
  PulseInteraction(const std::shared_ptr<const DotVector> &,
                   const std::shared_ptr<const Pulse>, const double, const double);
  virtual const ResultArray &evaluate(const int);

 private:
  std::shared_ptr<const Pulse> pulse;
  const double hbar;
  const double dt;
};

#endif
