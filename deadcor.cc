#include <cstdio>
#include <stdexcept>
#include <string>

#include "deadcor.hh"
#include "common.hh"

DeadCorTable::DeadCorTable(fitsfile *ff, int tm)
  : cache_idx(0)
{
  int status = 0;

  std::string hduname = std::string("DEADCOR") + std::to_string(tm);
  std::printf("  - Opening extension %s\n", hduname.c_str());
  move_fits_hdu(ff, hduname.c_str());

  long nrows;
  fits_get_num_rows(ff, &nrows, &status);
  check_fitsio_status(status);

  num = nrows;
  time.resize(num);
  deadc.resize(num);

  read_fits_column(ff, "TIME", TDOUBLE, nrows, &time[0]);
  read_fits_column(ff, "DEADC", TFLOAT, nrows, &deadc[0]);

  std::printf("    - successfully read %ld entries\n", nrows);
}

float DeadCorTable::interpolate(double t)
{
  // simple linear scan - should be ok if we're processing in time order
  while( time[cache_idx+1]>t and cache_idx>0 )
    --cache_idx;

  while( t>time[cache_idx+1] and cache_idx<int(time.size())-1 )
    ++cache_idx;

  if( t<time[cache_idx] || t>time[cache_idx+1] )
    throw std::runtime_error("Interpolated DEADCOR outside of time");

  float f = float((t - time[cache_idx]) /
                  (time[cache_idx+1] - time[cache_idx]));

  return deadc[cache_idx]*(1-f) + deadc[cache_idx+1]*f;
}
