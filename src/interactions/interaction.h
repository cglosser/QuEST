#ifndef INTERACTION_H
#define INTERACTION_H

#include <vector>

#include "../quantum_dot.h"

class Interaction {
 public:
  Interaction(const std::shared_ptr<const DotVector> &dots) : dots(dots){};
  double result(const int i) { return results[i]; }
  virtual void evaluate(const int);

 protected:
  std::shared_ptr<const DotVector> dots;
  std::vector<double> results;
};

#endif
