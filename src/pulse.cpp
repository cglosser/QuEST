#include "pulse.h"

Pulse::Pulse(const double amplitude, const double delay, const double width,
             const double freq, const Eigen::Vector3d &wavevector,
             const Eigen::Vector3d &polarization,
             const Configuration::REFERENCE_FRAME ref_frame) //add a fixed vs. rot argument
    : amplitude(amplitude),
      delay(delay),
      width(width),
      freq(freq),
      wavevector(wavevector),
      polarization(polarization.normalized()),
      ref_frame(ref_frame) //add a fixed vs rot argument
{
}

Eigen::Vector3d Pulse::operator()(const Eigen::Vector3d &r,
                                  const double t) const
{
  const double arg = wavevector.dot(r) - freq * (t - delay);

  if(ref_frame == Configuration::REFERENCE_FRAME::ROTATING) {
    return (amplitude / 2 * polarization) * gaussian(arg / width);  // * cos(arg); write an if statement that determines whether the cos(arg) term is included

  } else {
    return (amplitude / 2 * polarization) * gaussian(arg / width) * cos(arg);
  }

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

void set_reference_frame(Pulse &p, Configuration::REFERENCE_FRAME ref_frame)
{
  p.ref_frame = ref_frame;
}

Pulse read_pulse_config(const std::string &filename, const Configuration::REFERENCE_FRAME ref_frame) //read in config.ref_frame
{
  std::ifstream ifs(filename);
  if(!ifs) throw std::runtime_error("Could not open " + filename);

  Pulse p;
  ifs >> p; //find a way to include config.ref_frame in p
  set_reference_frame(p,ref_frame);

  return p;
}
