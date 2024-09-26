#ifndef DEADCOR_H
#define DEADCOR_H

#include <vector>
#include <fitsio.h>

class DeadCorTable
{
  public:
  DeadCorTable(fitsfile *ff, int tm);

  float interpolate(double t);

  size_t num;
  std::vector<double> time;
  std::vector<float> deadc;

  int cache_idx;
};

#endif
