#ifndef MODE_HH
#define MODE_HH

#include "geom.hh"

class ProjMode
{
public:
  // use source when it's in the current position?
  virtual bool sourceValid(Point ccdpt) const = 0;

  // get matrix to rotate photons or exposure
  virtual RotationMatrix rotationMatrix(double roll, Point delccd) const;

  // origin to use given source
  virtual Point origin(Point ccdpt) const;
};

// average field of view, in CCD coordinates
class ProjModeAverageFoV : public ProjMode
{
public:
  bool sourceValid(Point ccdpt) const;
};

// average FoV, in sky-relative coordinates
class ProjModeAverageFoVSky : public ProjMode
{
public:
  bool sourceValid(Point ccdpt) const;
  RotationMatrix rotationMatrix(double roll, Point delccd) const;
};

// mode for testing - just Det coordinates without tracking source
class ProjModeDet : public ProjMode
{
public:
  bool sourceValid(Point ccdpt) const;
  Point origin(Point ccdpt) const;
};


#endif
