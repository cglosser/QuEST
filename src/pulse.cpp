#include "pulse.h"

Pulse::Pulse(const double amplitude, const double delay, const double width,
             const double freq, const Eigen::Vector3d &wavevector,
             const Eigen::Vector3d &polarization)
    : amplitude(amplitude),
      delay(delay),
      width(width),
      freq(freq),
      wavevector(wavevector),
      polarization(polarization.normalized())
{
}

Eigen::Vector3d Pulse::operator()(const Eigen::Vector3d &r,
                                  const double t) const
{
  const double arg = wavevector.dot(r) - freq * (t - delay);
  return (amplitude / 2 * polarization) * gaussian(arg / width);  // * cos(arg);
}

std::ostream &operator<<(std::ostream &os, const Pulse &p)
{
  os << p.amplitude << " " << p.delay << " " << p.width << " " << p.freq << " "
     << p.wavevector.transpose() << " " << p.polarization.transpose();
  return os;
}

std::istream &operator>>(std::istream &is, Pulse &p)
{
  is >> p.amplitude >> p.delay >> p.width >> p.freq >> p.wavevector[0] >>
      p.wavevector[1] >> p.wavevector[2] >> p.polarization[0] >>
      p.polarization[1] >> p.polarization[2];
  return is;
}

Pulse read_pulse_config(const std::string &fname)
{
  std::ifstream ifs(fname);
  if(!ifs) throw std::runtime_error("Could not open " + fname);

  Pulse p;
  ifs >> p;

  return p;
}
