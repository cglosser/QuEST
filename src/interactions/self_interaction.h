#ifndef SELF_INTERACTION
#define SELF_INTERACTION

#include <Eigen/Dense>
#include "../integrator/history.h"
#include "interaction.h"

class SelfInteraction : public Interaction {
 public:
  SelfInteraction(const std::shared_ptr<const DotVector> &,
                  const std::shared_ptr<const Integrator::History<soltype>> &);
  virtual const ResultArray &evaluate(const int);

 private:
  std::shared_ptr<const Integrator::History<soltype>> history;
};

#endif
