#ifndef MASK_HH
#define MASK_HH

#include <string>
#include <vector>

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
  Mask(const std::string& filename);

  CoordVecVec maskcoords;
};

#endif
