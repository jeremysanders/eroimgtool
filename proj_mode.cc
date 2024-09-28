#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "common.hh"
#include "proj_mode.hh"

// FIXME: assumption centre of coords is 192,192

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

void ProjModeAverageFoV::message() const
{
  std::printf("Projection mode\n");
  std::printf("  - fov: source-relative detector coordinates for std FoV\n");
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

void ProjModeAverageFoVSky::message() const
{
  std::printf("Projection mode\n");
  std::printf("  - fov_sky: source-relative sky-rotated detector coordinates for std FoV\n");
}

////////////////////////////////////////////////////////////////////


Point ProjModeDet::origin(Point ccdpt) const
{
  return Point(192,192);
}

bool ProjModeDet::sourceValid(Point ccdpt) const
{
  return true;
}

void ProjModeDet::message() const
{
  std::printf("Projection mode\n");
  std::printf("  - det: non-relative detector coordinates\n");
}

////////////////////////////////////////////////////////////////////

ProjModeRadial::ProjModeRadial(const std::vector<float>& args)
{
  if(args.size() != 2)
    throw std::runtime_error("Two parameters are required for radial projection (rin,rout)");

  rin = args[0];
  rout = args[1];
}

bool ProjModeRadial::sourceValid(Point ccdpt) const
{
  float rad = std::sqrt(sqr(ccdpt.x-192) + sqr(ccdpt.y-192));
  return (rad >= rin) && (rad < rout);
}

void ProjModeRadial::message() const
{
  std::printf("Projection mode\n");
  std::printf("  - radial: radial range of detector (%g to %g pix)\n",
              rin, rout);
}

////////////////////////////////////////////////////////////////////

void ProjModeRadialSym::message() const
{
  std::printf("Projection mode\n");
  std::printf("  - radial symmetric: radial range of detector (%g to %g pix)\n",
              rin, rout);
}

RotationMatrix ProjModeRadialSym::rotationMatrix(double roll,
                                                 Point delccd) const
{
  float theta = -std::atan2(delccd.y, delccd.x);
  float c = std::cos(theta);
  float s = std::sin(theta);
  return RotationMatrix(c, -s, s, c);
}
