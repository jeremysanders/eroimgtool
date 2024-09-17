#ifndef BADPIX_HH
#define BADPIX_HH

#include <vector>
#include <fitsio.h>
#include "image.hh"
#include "geom.hh"

class BadPixTable
{
public:
  BadPixTable(fitsfile *ff, int tm);

  const PolyVec& getPolyMask(double t);

private:
  void buildMask(double t);

private:
  size_t num_entries;
  std::vector<int> rawx, rawy, yextent;
  std::vector<double> timemin, timemax;

  // times where bad pixel table changes
  std::vector<double> tedge;

  int cache_ti;
  PolyVec cache_poly;
  Image<int> cache_mask;
};

#endif
