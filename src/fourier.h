#ifndef FOURIER_H
#define FOURIER_H

#include <fftw3.h>

struct TransformPair {
  fftw_plan forward, backward;
  ~TransformPair()
  {
    if(forward) fftw_destroy_plan(forward);
    if(backward) fftw_destroy_plan(backward);
  }
};

#endif
