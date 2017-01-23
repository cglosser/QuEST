#ifndef INTERACTION_H
#define INTERACTION_H

#include <Eigen/Dense>

#include "../quantum_dot.h"

class Interaction {
 public:
  typedef Eigen::Array<cmplx, Eigen::Dynamic, 1> ResultArray;

  Interaction(const std::shared_ptr<const DotVector> &dots)
      : dots(dots), results(dots->size()){};
  cmplx result(const int i) { return results[i]; }
  virtual void evaluate(const int) = 0;

 protected:
  std::shared_ptr<const DotVector> dots;
  ResultArray results;
};

#endif
