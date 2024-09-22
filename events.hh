#ifndef EVENTS_HH
#define EVENTS_HH

#include <vector>
#include <fitsio.h>

class Events
{
public:
  Events(fitsfile *ff);
  void filter_tm(int tm);
  void filter_pi(float pimin, float pimax);

private:
  void do_filter(const std::vector<size_t>& sel);

public:
  size_t num_entries;
  std::vector<short> rawx, rawy, tm_nr;
  std::vector<double> ra, dec, time;
  std::vector<float> pi, subx, suby;

  // combined rawx+subx, rawy+suby
  std::vector<float> ccdx, ccdy;
};

#endif
