#include <cstdio>
#include <cstring>
#include <string>
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
