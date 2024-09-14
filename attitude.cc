#include <cstdio>
#include "attitude.hh"
#include "common.hh"

AttitudeTable::AttitudeTable(fitsfile *ff, int tm)
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
}
