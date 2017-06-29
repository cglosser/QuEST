#ifndef INTERACTION_H
#define INTERACTION_H

#include <memory>
#include <Eigen/Dense>

#include "../magnetic_particle.h"

class Interaction {
 public:
  typedef std::vector<Eigen::Vector3d> ResultArray;

  Interaction(const std::shared_ptr<const DotVector> &dots)
      : dots(dots), results(dots->size()){};
  const Eigen::Vector3d &operator[](const int i) const { return results[i]; }
  virtual const ResultArray &evaluate(const int) = 0;

 protected:
  std::shared_ptr<const DotVector> dots;
  ResultArray results;
};

#endif
