#include <cstdio>
#include "badpix.hh"
#include "common.hh"

BadPixTable::BadPixTable(fitsfile *ff, int tm)
{
  int status = 0;

  char hduname[16];
  std::sprintf(hduname, "BADPIX%d", tm);

  std::printf("  Opening bad pixel extension %s\n", hduname);
  fits_movnam_hdu(ff, BINARY_TBL, hduname, 0, &status);
  check_fitsio_status(status);

  // get indices of columns
  int c_rawx, c_rawy, c_yextent, c_timemin, c_timemax;
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("RAWX"),
                  &c_rawx, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("RAWY"),
                  &c_rawy, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("YEXTENT"),
                  &c_yextent, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("TIMEMIN"),
                  &c_timemin, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("TIMEMAX"),
                  &c_timemax, &status);

  long nrows;
  fits_get_num_rows(ff, &nrows, &status);

  check_fitsio_status(status);

  num = nrows;
  rawx.resize(num);
  rawy.resize(num);
  yextent.resize(num);
  timemin.resize(num);
  timemax.resize(num);

  // read in table
  fits_read_col(ff, TINT, c_rawx, 1, 1, nrows, 0,
                &rawx[0], 0, &status);
  fits_read_col(ff, TINT, c_rawy, 1, 1, nrows, 0,
                &rawy[0], 0, &status);
  fits_read_col(ff, TINT, c_yextent, 1, 1, nrows, 0,
                &yextent[0], 0, &status);
  fits_read_col(ff, TDOUBLE, c_timemin, 1, 1, nrows, 0,
                &timemin[0], 0, &status);
  fits_read_col(ff, TDOUBLE, c_timemax, 1, 1, nrows, 0,
                &timemax[0], 0, &status);
  check_fitsio_status(status);

  std::printf("    - successfully read %ld entries\n", nrows);
}
