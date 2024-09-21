#ifndef MODE_HH
#define MODE_HH

#include <array>

#include "geom.hh"

class Mode
{
public:
  // use source when it's in the current position?
  virtual bool source_valid(Point ccdpt) const = 0;

  // get matrix to rotate photons or exposure
  virtual std::array<float,4> rotation_matrix(double roll, Point delccd) const;

  // origin to use given source
  virtual Point origin(Point ccdpt) const;
};

// average field of view, in CCD coordinates
class ModeAverageFoV : public Mode
{
public:
  bool source_valid(Point ccdpt) const;
};

// average FoV, in sky-relative coordinates
class ModeAverageFoVSky : public Mode
{
public:
  bool source_valid(Point ccdpt) const;
  std::array<float,4> rotation_matrix(double roll, Point delccd) const;
};


#endif
