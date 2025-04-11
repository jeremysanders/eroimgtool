#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "common.hh"
#include "proj_mode.hh"

// FIXME: assumption centre of coords is 192,192

Matrix2 ProjMode::rotationMatrix(double roll, Point delccd) const
{
  return Matrix2();
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
  std::printf("  - fov: source-relative det coords for source in std FoV\n");
}

////////////////////////////////////////////////////////////////////

bool ProjModeAverageFoVSky::sourceValid(Point ccdpt) const
{
  float rad2 = sqr(ccdpt.x-192) + sqr(ccdpt.y-192);

  // FIXME: check radius
  return rad2 < sqr(192);
}

Matrix2 ProjModeAverageFoVSky::rotationMatrix(double roll, Point delccd) const
{
  float c = std::cos((270-roll)*DEG2RAD);
  float s = std::sin((270-roll)*DEG2RAD);
  return Matrix2(c, -s, s, c);
}

void ProjModeAverageFoVSky::message() const
{
  std::printf("Projection mode\n");
  std::printf("  - fov_sky: source-relative sky-rotated det coords for source in std FoV\n");
}

////////////////////////////////////////////////////////////////////

bool ProjModeAverageFull::sourceValid(Point ccdpt) const
{
  return true;
}

void ProjModeAverageFull::message() const
{
  std::printf("Projection mode\n");
  std::printf("  - full: source-relative detector coordinates\n");
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
  if(args.size() != 2 && args.size() != 4)
    throw std::runtime_error("2 or 4 parameters are required for radial projection (rin,rout) or (rin,rout,cx,cy)");

  rin = args[0];
  rout = args[1];

  if(args.size() == 4)
    {
      cx = args[2];
      cy = args[3];
    }
  else
    {
      cx = 192;
      cy = 192;
    }
}

bool ProjModeRadial::sourceValid(Point ccdpt) const
{
  float rad = std::sqrt(sqr(ccdpt.x-cx) + sqr(ccdpt.y-cy));
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

Matrix2 ProjModeRadialSym::rotationMatrix(double roll, Point delccd) const
{
  float theta = -std::atan2(delccd.y, delccd.x);
  float s = std::sin(theta);
  float c = std::cos(theta);
  return Matrix2(c, -s, s, c);
}

////////////////////////////////////////////////////////////////////

ProjModeBox::ProjModeBox(const std::vector<float>& args)
{
  if(args.size() != 4)
    throw std::runtime_error("Four parameters are required for box (x1,y1,x2,y2)");

  x1 = args[0];
  y1 = args[1];
  x2 = args[2];
  y2 = args[3];
}

bool ProjModeBox::sourceValid(Point ccdpt) const
{
  return (ccdpt.x >= x1 && ccdpt.y >= y1 &&
          ccdpt.x <  x2 && ccdpt.y <  y2);
}

void ProjModeBox::message() const
{
  std::printf("Projection mode\n");
  std::printf("  - box: range (x=%g:%g,y=%g:%g pix)\n",
              x1, x2, y1, y2);
}
