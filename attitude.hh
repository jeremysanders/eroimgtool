#ifndef ATTITUDE_HH
#define ATTITUDE_HH

#include <vector>
#include <fitsio.h>

class AttitudeTable
{
  public:
  AttitudeTable(fitsfile *ff, int tm);

  size_t num;
  std::vector<double> time, ra, dec, roll;
};

#endif
