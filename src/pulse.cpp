#include "pulse.h"

Pulse::Pulse(const double a, const double t0, const double sigma,
             const double freq, const Eigen::Vector3d &k,
             const Eigen::Vector3d &e0)
    : amplitude(a),
      delay(t0),
      width(sigma),
      frequency(freq),
      wavevector(k),
      polarization(e0.normalized())
{
}

Eigen::Vector3d Pulse::operator()(const Eigen::Vector3d &r, const double t) const
{
  const double arg = wavevector.dot(r) - frequency * (t - delay);
  return (amplitude * polarization) * gaussian(arg / width) * cos(arg);
}

std::ostream &operator<<(std::ostream &os, const Pulse &p)
{
  os << p.amplitude << " " << p.delay << " " << p.width << " " << p.frequency
     << " " << p.wavevector.transpose() << " " << p.polarization.transpose();
  return os;
}

std::istream &operator>>(std::istream &is, Pulse &p)
{
  is >> p.amplitude >> p.delay >> p.width >> p.frequency >> p.wavevector[0] >>
      p.wavevector[1] >> p.wavevector[2] >> p.polarization[0] >>
      p.polarization[1] >> p.polarization[2];
  return is;
}
