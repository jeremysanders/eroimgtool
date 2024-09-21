#include <cstdio>
#include <string>
#include "gti.hh"
#include "common.hh"

GTITable::GTITable(fitsfile *ff, int tm)
{
  int status = 0;

  std::string hduname = std::string("GTI") + std::to_string(tm);
  std::printf("  - Opening GTI extension %s\n", hduname.c_str());
  move_fits_hdu(ff, hduname.c_str());

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
