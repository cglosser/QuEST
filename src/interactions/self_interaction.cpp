#include "self_interaction.h"

SelfInteraction::SelfInteraction(
    const std::shared_ptr<const DotVector> &dots,
    const std::shared_ptr<const Integrator::History<soltype>> &history)
    : Interaction(dots), history(history)
{
}

const Interaction::ResultArray &SelfInteraction::evaluate(const int time_idx)
{
  // result = history[particle][dt][0] * 1/3;
  for(int i = 0; i < dots->size(); ++i) {
    results[i] = history->array[i][time_idx][0] / 3;
    return results;
  }
}
