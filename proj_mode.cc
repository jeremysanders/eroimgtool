#include <cmath>
#include "common.hh"
#include "proj_mode.hh"

RotationMatrix ProjMode::rotationMatrix(double roll, Point delccd) const
{
  return RotationMatrix();
}

Point ProjMode::origin(Point ccdpt) const
{
  return ccdpt;
}

////////////////////////////////////////////////////////////////////

bool ProjModeAverageFoV::sourceValid(Point ccdpt) const
{
  float rad2 = sqr(ccdpt.x-192) + sqr(ccdpt.y-192);

  // FIXME: check radius
  return rad2 < sqr(192);
}

////////////////////////////////////////////////////////////////////

bool ProjModeAverageFoVSky::sourceValid(Point ccdpt) const
{
  float rad2 = sqr(ccdpt.x-192) + sqr(ccdpt.y-192);

  // FIXME: check radius
  return rad2 < sqr(192);
}

RotationMatrix ProjModeAverageFoVSky::rotationMatrix(double roll,
                                                     Point delccd) const
{
  float c = std::cos((270-roll)*DEG2RAD);
  float s = std::sin((270-roll)*DEG2RAD);
  return RotationMatrix(c, -s, s, c);
}



