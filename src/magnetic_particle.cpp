#include "magnetic_particle.h"

MagneticParticle::MagneticParticle(const Eigen::Vector3d &pos,
                                   const double alpha, const double gamma0,
                                   const double sat_mag, const soltype &mag)
    : pos(pos),
      alpha(alpha),
      gamma0(gamma0), 
      mag(mag),
      sat_mag(sat_mag)
{
}

soltype MagneticParticle::llg_rhs(const soltype &mag,
                                  const Eigen::Vector3d &hfield)
{
  const Eigen::Vector3d mxh = mag.cross(hfield);
  const double gamma = gamma0 / ( 1 + std::pow(alpha,2));

  return -gamma * mxh - gamma * alpha / sat_mag * mag.cross(mxh);
}

const Eigen::Vector3d separation(const MagneticParticle &mp1,
                                 const MagneticParticle &mp2)
{
  return mp2.pos - mp1.pos;
}

rhs_func_vector rhs_functions(const DotVector &dots)
{
  rhs_func_vector funcs(dots.size());

  using std::placeholders::_1;
  using std::placeholders::_2;
  std::transform(dots.begin(), dots.end(), funcs.begin(),
                 [](const MagneticParticle &mp) {
                   return std::bind(&MagneticParticle::llg_rhs, mp, _1, _2);
                 });
  return funcs;
}

std::ostream &operator<<(std::ostream &os, const MagneticParticle &mp)
{
  os << mp.pos.transpose() << " " << mp.alpha << " " << mp.gamma0 << " "
     << mp.sat_mag << " " << mp.mag.transpose() << std::endl;
  return os;
}

std::istream &operator>>(std::istream &is, MagneticParticle &mp)
{
  is >> mp.pos[0] >> mp.pos[1] >> mp.pos[2] >> mp.alpha >> mp.gamma0 >>
      mp.sat_mag >> mp.mag[0] >> mp.mag[1] >> mp.mag[2];
  return is;
}

DotVector import_dots(const std::string &fname)
{
  std::ifstream ifs(fname);
  if(!ifs) throw std::runtime_error("Could not open " + fname);

  std::istream_iterator<MagneticParticle> in_iter(ifs), eof;
  return DotVector(in_iter, eof);
}
