#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "attitude.hh"
#include "common.hh"

AttitudeTable::AttitudeTable(fitsfile *ff, int tm)
  : cache_idx(0)
{
  int status = 0;

  std::string hduname = std::string("CORRATT") + std::to_string(tm);
  std::printf("  - Opening attitude extension %s\n", hduname.c_str());
  move_fits_hdu(ff, hduname.c_str());

  long nrows;
  fits_get_num_rows(ff, &nrows, &status);
  check_fitsio_status(status);

  num = nrows;
  time.resize(nrows);
  ra.resize(num);
  dec.resize(num);
  roll.resize(num);

  read_fits_column(ff, "TIME", TDOUBLE, nrows, &time[0]);
  read_fits_column(ff, "RA", TDOUBLE, nrows, &ra[0]);
  read_fits_column(ff, "DEC", TDOUBLE, nrows, &dec[0]);
  read_fits_column(ff, "ROLL", TDOUBLE, nrows, &roll[0]);

  std::printf("    - successfully read %ld entries\n", nrows);
  std::printf("    - time from %.0f to %.0f\n", time[0], time[nrows-1]);
}

std::tuple<double, double, double> AttitudeTable::interpolate(double t)
{
  while( time[cache_idx+1]>t and cache_idx>0 )
    --cache_idx;

  while( t>time[cache_idx+1] and cache_idx<int(time.size())-1 )
    ++cache_idx;

  if( t<time[cache_idx] || t>time[cache_idx+1] )
    throw std::runtime_error("Interpolated attitude outside of time");

  double f = (t - time[cache_idx]) / (time[cache_idx+1] - time[cache_idx]);

  double ra_int = ra[cache_idx]*(1-f) + ra[cache_idx+1]*f;
  double dec_int = dec[cache_idx]*(1-f) + dec[cache_idx+1]*f;

  // interpolate sin and cos, and use atan2 to get back angle
  double s1 = std::sin(roll[cache_idx]*DEG2RAD);
  double c1 = std::cos(roll[cache_idx]*DEG2RAD);
  double s2 = std::sin(roll[cache_idx+1]*DEG2RAD);
  double c2 = std::cos(roll[cache_idx+1]*DEG2RAD);
  double roll_int = std::atan2(s1*(1-f)+s2*f, c1*(1-f)+c2*f)*RAD2DEG;

  // combine to return
  return std::make_tuple(ra_int, dec_int, roll_int);
}
