#include "bloch_system.h"

BlochSystem::BlochSystem(const Weights pc,
                         std::vector<QuantumDot> ds, const int iorder,
                         const size_t nsteps)
    : num_steps(nsteps + 1),
      num_dots(static_cast<int>(ds.size())),
      integrator(pc),
      history(boost::extents[ds.size()]
                            [HistoryArray::extent_range(-pc.width(), num_steps)]
                            [2]),
      interactions(iorder, ds),
      rabi_freqs(ds.size()),
      now(1),
      dt(1)
{
  dots.swap(ds);

  for(int dot_idx = 0; dot_idx < static_cast<int>(dots.size()); ++dot_idx) {
    for(int i = -pc.width(); i <= 0; ++i) {
      history[dot_idx][i][0] = matrix_elements(1, 0);
      history[dot_idx][i][1] = matrix_elements(0, 0);
    }
  }
}

void BlochSystem::step()
{
  assert(now < num_steps);

  const int start = now - integrator.width();
  const double time = now * dt, norm = 0.00489575834889;
  std::vector<double> rabi(
      num_dots, norm * gaussian((time - 1024) / 256) * cos(0.05 * time));

  for(int dot_idx = 0; dot_idx < num_dots; ++dot_idx) {
    // Predictor
    for(int h = 0; h < integrator.width(); ++h) {
      history[dot_idx][now][0] +=
          history[dot_idx][start + h][0] * integrator.predictor_coefs(h, 0) +
          history[dot_idx][start + h][1] * integrator.predictor_coefs(h, 1) *
              dt;
    }

    // Estimator
    convolve_currents(now);
    history[dot_idx][now][1] = dots[dot_idx].liouville_rhs(
        history[dot_idx][now][0], rabi_freqs[dot_idx]);

    for(int m = 0; m < 10; ++m) {
      // Corrector
      history[dot_idx][now][0] =
          integrator.future_coef * history[dot_idx][now][1] * dt;

      for(int h = 0; h < integrator.width(); ++h) {
        history[dot_idx][now][0] +=
            history[dot_idx][start + h][0] * integrator.corrector_coefs(h, 0) +
            history[dot_idx][start + h][1] * integrator.corrector_coefs(h, 1) *
                dt;
      }

      // Estimator
      convolve_currents(now);
      history[dot_idx][now][1] = dots[dot_idx].liouville_rhs(
          history[dot_idx][now][0], rabi_freqs[dot_idx]);
    }
  }

  now++;
}

void BlochSystem::convolve_currents(const int time_idx)
{
  std::fill(rabi_freqs.begin(), rabi_freqs.end(), 0);

  for(size_t src = 0; src < dots.size() - 1; ++src) {
    for(size_t obs = src + 1; obs < dots.size(); ++obs) {
      const size_t idx = interactions.coord2idx(src, obs);
      const int s = time_idx - interactions.floor_delays[idx] - 1;

      for(int i = 0; i <= interactions.interp_order; ++i) {
        if(s - i < history.index_bases()[1]) continue;
        rabi_freqs[src] += polarization(history[obs][s - i][0]) *
                           interactions.coefficients[idx][i];
        rabi_freqs[obs] += polarization(history[src][s - i][0]) *
                           interactions.coefficients[idx][i];
      }
    }
  }
}
