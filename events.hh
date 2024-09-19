#ifndef EVENTS_HH
#define EVENTS_HH

#include <vector>
#include <fitsio.h>

class Events
{
public:
  Events(fitsfile *ff, int tm);

public:
  size_t num_entries;
  std::vector<short> rawx, rawy;
  std::vector<double> ra, dec, time;
  std::vector<float> pi, subx, suby;

  // combined rawx+subx, rawy+suby
  std::vector<float> ccdx, ccdy;
};

#endif
