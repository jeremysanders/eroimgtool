#include <cmath>
#include "common.hh"
#include "mode.hh"

std::array<float,4> Mode::rotation_matrix(double roll, Point delccd) const
{
  return {1,0,0,1};
}

Point Mode::origin(Point ccdpt) const
{
  return ccdpt;
}

////////////////////////////////////////////////////////////////////

bool ModeAverageFoV::source_valid(Point ccdpt) const
{
  float rad2 = sqr(ccdpt.x-192) + sqr(ccdpt.y-192);

  // FIXME: check radius
  return rad2 < sqr(192);
}

////////////////////////////////////////////////////////////////////

bool ModeAverageFoVSky::source_valid(Point ccdpt) const
{
  float rad2 = sqr(ccdpt.x-192) + sqr(ccdpt.y-192);

  // FIXME: check radius
  return rad2 < sqr(192);
}

std::array<float,4> ModeAverageFoVSky::rotation_matrix(double roll,
                                                       Point delccd) const
{
  float c = std::cos((270-roll)*DEG2RAD);
  float s = std::sin((270-roll)*DEG2RAD);
  return {c,-s,s,c};
}



