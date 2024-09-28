#include <algorithm>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "gti.hh"
#include "common.hh"

GTITable::GTITable(fitsfile *ff, int tm)
{
  int status = 0;

  // try GTIx and then STDGTI
  char hdu[32];
  std::sprintf(hdu, "GTI%d", tm);
  fits_movnam_hdu(ff, ANY_HDU, hdu, 0, &status);
  if(status != 0)
    {
      // try STDGTI hdu
      status = 0;
      fits_clear_errmsg();
      std::strcpy(hdu, "STDGTI");
      fits_movnam_hdu(ff, ANY_HDU, hdu, 0, &status);
    }
  check_fitsio_status(status);

  std::printf("  - Opening GTI extension %s\n", hdu);
  long nrows;
  fits_get_num_rows(ff, &nrows, &status);
  check_fitsio_status(status);

  num = nrows;
  start.resize(num);
  stop.resize(num);

  read_fits_column(ff, "START", TDOUBLE, nrows, &start[0]);
  read_fits_column(ff, "STOP", TDOUBLE, nrows, &stop[0]);

  std::printf("    - successfully read %ld entries\n", nrows);
}

void GTITable::operator&=(const GTITable& o)
{
  struct timesign
  {
    double t;
    int sign;
    bool operator<(const timesign& o) const { return t<o.t; }
  };

  // collect together time periods
  std::vector<timesign> times;
  for(auto t : start)
    times.emplace_back(timesign{t,+1});
  for(auto t : stop)
    times.emplace_back(timesign{t,-1});
  for(auto t : o.start)
    times.emplace_back(timesign{t,+1});
  for(auto t : o.stop)
    times.emplace_back(timesign{t,-1});

  // sort them in time order
  std::sort(times.begin(), times.end());

  start.clear();
  stop.clear();
  int ct = 0;
  for(auto& ts : times)
    {
      // keep track of "stack" of times, and record when two are active
      ct += ts.sign;
      if(ct == 2 && ts.sign==1)
        {
          if(!stop.empty() && stop.back() == ts.t)
            stop.pop_back();
          else
            start.emplace_back(ts.t);
        }
      else if(ct == 1 && ts.sign==-1)
        {
          if(!start.empty() && start.back() == ts.t)
            start.pop_back();
          else
            stop.emplace_back(ts.t);
        }
    }
  num = start.size();

  if(start.size() != stop.size())
    throw std::runtime_error("internal error merging GTIs");
}
