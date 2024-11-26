#ifndef DETMAP_HH
#define DETMAP_HH

#include <vector>
#include <fitsio.h>
#include "image.hh"

class DetMap
{
public:
  DetMap(fitsfile *ff, int tm, bool detmapmask);

  const Image<float>& getMap(double t) { checkCache(t); return cache_map; }

private:
  void checkCache(double t);
  void buildMapImage(double t);
  void readDetmapMask(int tm);

private:
  size_t num_entries;
  std::vector<int> rawx, rawy, yextent;
  std::vector<double> timemin, timemax;

  // times where bad pixel table changes
  std::vector<double> tedge;

  int cache_ti;
  Image<float> init_map, cache_map;
};

#endif
