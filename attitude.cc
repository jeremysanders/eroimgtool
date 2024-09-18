#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "attitude.hh"
#include "common.hh"

AttitudeTable::AttitudeTable(fitsfile *ff, int tm)
  : cache_idx(0)
{
  int status = 0;

  char hduname[16];
  std::sprintf(hduname, "CORRATT%d", tm);

  std::printf("  Opening attitude extension %s\n", hduname);
  fits_movnam_hdu(ff, BINARY_TBL, hduname, 0, &status);
  check_fitsio_status(status);

  // get indices of columns
  int c_time, c_ra, c_dec, c_roll;
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("TIME"),
                  &c_time, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("RA"),
                  &c_ra, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("DEC"),
                  &c_dec, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("ROLL"),
                  &c_roll, &status);

  long nrows;
  fits_get_num_rows(ff, &nrows, &status);

  check_fitsio_status(status);

  num = nrows;
  time.resize(nrows);
  ra.resize(num);
  dec.resize(num);
  roll.resize(num);

  // read in table
  fits_read_col(ff, TDOUBLE, c_time, 1, 1, nrows, 0,
                &time[0], 0, &status);
  fits_read_col(ff, TDOUBLE, c_ra, 1, 1, nrows, 0,
                &ra[0], 0, &status);
  fits_read_col(ff, TDOUBLE, c_dec, 1, 1, nrows, 0,
                &dec[0], 0, &status);
  fits_read_col(ff, TDOUBLE, c_roll, 1, 1, nrows, 0,
                &roll[0], 0, &status);
  check_fitsio_status(status);

  std::printf("    - successfully read %ld entries\n", nrows);
  std::printf("    - time from %.0f to %.0f\n", time[0], time[nrows-1]);
}

std::tuple<double, double, double> AttitudeTable::interpolate(double t)
{
  while( time[cache_idx+1]>t and cache_idx>0 )
    --cache_idx;
  while( t>time[cache_idx] and cache_idx<int(time.size())-1 )
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
  double roll_int = std::atan2(s1*(1-f)+s2*f, c1*(1-f)+c2*f);

  // combine to return
  return std::make_tuple(ra_int, dec_int, roll_int);
}
