#ifndef GTI_HH
#define GTI_HH

#include <vector>
#include <fitsio.h>

class GTITable
{
  public:
  GTITable(fitsfile *ff, int tm);

  size_t num;
  std::vector<double> start, stop;
};

#endif
