#ifndef DETMAP_HH
#define DETMAP_HH

#include <string>
#include <vector>
#include <fitsio.h>
#include "image.hh"

class DetMap
{
public:
  // create detector mask
  // tm: TM number
  // detmapmask: mask using standard CALDB detector mask
  // shadowmask: mask out bottom pixels in readout shadow area
  DetMap(int _tm, bool detmapmask, bool shadowmask);

  // read table from open fits file
  void read(fitsfile *ff);

  // read table from fits file given
  void read(const std::string& filename);

  const Image<float>& getMap(double t) { checkCache(t); return cache_map; }

private:
  void checkCache(double t);
  void buildMapImage(double t);
  void readDetmapMask(int tm);

private:
  int tm;

  size_t num_entries;
  std::vector<int> rawx, rawy, yextent;
  std::vector<double> timemin, timemax;

  // times where bad pixel table changes
  std::vector<double> tedge;

  int cache_ti;
  Image<float> init_map, cache_map;
};

#endif
