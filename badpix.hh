#ifndef BADPIX_HH
#define BADPIX_HH

#include <vector>
#include <fitsio.h>

class BadPixTable
{
  public:
  BadPixTable(fitsfile *ff, int tm);

  size_t num;
  std::vector<int> rawx, rawy, yextent;
  std::vector<double> timemin, timemax;
};

#endif
