#ifndef ATTITUDE_HH
#define ATTITUDE_HH

#include <tuple>
#include <vector>

#include <fitsio.h>

class AttitudeTable
{
  public:
  AttitudeTable(fitsfile *ff, int tm);

  // return ra, dec, roll, interpolated for a time given
  std::tuple<double, double, double> interpolate(double t);

  size_t num;
  std::vector<double> time, ra, dec, roll;

  int cache_idx;
};

#endif
