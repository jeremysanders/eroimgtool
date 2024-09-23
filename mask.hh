#ifndef MASK_HH
#define MASK_HH

#include <string>
#include <vector>

#include "coords.hh"
#include "geom.hh"

struct Coord
{
  Coord() {}
  Coord(double _lon, double _lat) : lon(_lon), lat(_lat) {}
  double lon, lat;
};

typedef std::vector<Coord> CoordVec;
typedef std::vector<CoordVec> CoordVecVec;

class Mask
{
public:
  Mask();
  Mask(const std::string& filename);

  PolyVec as_ccd_poly(const CoordConv& cc) const;

private:
  CoordVecVec maskcoords;
};

#endif
