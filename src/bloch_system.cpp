#include "bloch_system.h"

BlochSystem::BlochSystem(const PredictorCorrector pc,
                         std::vector<QuantumDot> ds, const size_t nsteps)
    : integrator(pc),
      history(
          boost::extents[ds.size()]
                        [HistoryArray::extent_range(-pc.num_steps(), nsteps)]
                        [2]),
      num_dots(static_cast<int>(ds.size())),
      num_steps(nsteps)
{
  dots.swap(ds);

  for(int dot_idx = 0; dot_idx < static_cast<int>(dots.size()); ++dot_idx) {
    for(int i = -pc.num_steps(); i < 0; ++i) {
      history[dot_idx][i][0] = matrix_elements(1, 0);
      history[dot_idx][i][1] = matrix_elements(0, 0);
    }
  }

}
