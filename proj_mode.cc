#include <cmath>
#include "common.hh"
#include "proj_mode.hh"

std::array<float,4> ProjMode::rotation_matrix(double roll, Point delccd) const
{
  return {1,0,0,1};
}

Point ProjMode::origin(Point ccdpt) const
{
  return ccdpt;
}

////////////////////////////////////////////////////////////////////

bool ProjModeAverageFoV::source_valid(Point ccdpt) const
{
  float rad2 = sqr(ccdpt.x-192) + sqr(ccdpt.y-192);

  // FIXME: check radius
  return rad2 < sqr(192);
}

////////////////////////////////////////////////////////////////////

bool ProjModeAverageFoVSky::source_valid(Point ccdpt) const
{
  float rad2 = sqr(ccdpt.x-192) + sqr(ccdpt.y-192);

  // FIXME: check radius
  return rad2 < sqr(192);
}

std::array<float,4> ProjModeAverageFoVSky::rotation_matrix(double roll,
                                                           Point delccd) const
{
  float c = std::cos((270-roll)*DEG2RAD);
  float s = std::sin((270-roll)*DEG2RAD);
  return {c,-s,s,c};
}



