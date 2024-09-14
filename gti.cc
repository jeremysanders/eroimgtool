#include <cstdio>
#include "gti.hh"
#include "common.hh"

GTITable::GTITable(fitsfile *ff, int tm)
{
  int status = 0;

  char hduname[16];
  std::sprintf(hduname, "GTI%d", tm);

  std::printf("  Opening GTI extension %s\n", hduname);
  fits_movnam_hdu(ff, BINARY_TBL, hduname, 0, &status);
  check_fitsio_status(status);

  // get indices of columns
  int c_start, c_stop;
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("START"),
                  &c_start, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("STOP"),
                  &c_stop, &status);

  long nrows;
  fits_get_num_rows(ff, &nrows, &status);

  check_fitsio_status(status);

  num = nrows;
  start.resize(num);
  stop.resize(num);

  // read in table
  fits_read_col(ff, TDOUBLE, c_start, 1, 1, nrows, 0,
                &start[0], 0, &status);
  fits_read_col(ff, TDOUBLE, c_stop, 1, 1, nrows, 0,
                &stop[0], 0, &status);
  check_fitsio_status(status);

  std::printf("    - successfully read %ld entries\n", nrows);
}
