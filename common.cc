#include <cstring>
#include <string>
#include <stdexcept>
#include <fitsio.h>

#include "common.hh"

void check_fitsio_status(int status)
{
  if(status != 0)
    {
      std::string comb;
      char msg[81];
      msg[0] = '\0';

      for(;;)
        {
          int ret = fits_read_errmsg(msg);
          if(ret==0) break;
          comb += std::string(msg) + '\n';
        }
      
      throw std::runtime_error(comb);
    }
}

void read_fits_column(fitsfile *ff, const char* name, int type, long nrows,
                      void *retn)
{
  char col[80];
  std::strcpy(col, name);

  int status = 0;
  int cidx;
  fits_get_colnum(ff, CASEINSEN, col, &cidx, &status);
  check_fitsio_status(status);
  fits_read_col(ff, type, cidx, 1, 1, nrows, 0, retn, 0, &status);
  check_fitsio_status(status);
}

void move_fits_hdu(fitsfile* ff, const char* name)
{
  char tname[80];
  std::strcpy(tname, name);
  int status = 0;
  fits_movnam_hdu(ff, ANY_HDU, tname, 0, &status);
  check_fitsio_status(status);
}
