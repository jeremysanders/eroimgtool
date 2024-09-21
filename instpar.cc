#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <stdexcept>

#include <fitsio.h>

#include "common.hh"
#include "instpar.hh"

std::string lookup_cal(const std::string& subdir, const std::string& cmpt)
{
  char* caldb = getenv("CALDB");
  if(caldb == nullptr)
    throw std::runtime_error("CALDB environment variable not set");

  std::string root = std::string(caldb) + "/data/erosita/" + subdir;
  std::string idx_fname = root + "/caldb.indx";

  int status = 0;
  fitsfile* ff;

  fits_open_file(&ff, idx_fname.c_str(), READONLY, &status);
  check_fitsio_status(status);

  fits_movnam_hdu(ff, ANY_HDU, const_cast<char*>("CIF"), 0, &status);
  check_fitsio_status(status);

  // get columns and number of rows
  int c_cal_file, c_cal_cnam, c_cal_qual;
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("CAL_FILE"),
                  &c_cal_file, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("CAL_CNAM"),
                  &c_cal_cnam, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("CAL_QUAL"),
                  &c_cal_qual, &status);
  long nrows;
  fits_get_num_rows(ff, &nrows, &status);
  check_fitsio_status(status);

  char tempstr[256];
  char* tempstrarr[1] = {tempstr};
  int qual;

  // find first matching row in caldb index
  std::string res;
  for(long row=1; row <= nrows; ++row)
    {
      fits_read_col(ff, TSTRING, c_cal_cnam, row, 1, 1, 0,
                    &tempstrarr, 0, &status);
      fits_read_col(ff, TINT, c_cal_qual, row, 1, 1, 0,
                    &qual, 0, &status);
      check_fitsio_status(status);
      std::string cnam(tempstr);
      if(cnam == cmpt && qual == 0)
        {
          fits_read_col(ff, TSTRING, c_cal_file, row, 1, 1, 0,
                        &tempstrarr, 0, &status);
          check_fitsio_status(status);
          res = tempstr;
          break;
        }
    }

  if(res == "")
    throw std::runtime_error("Could not find calibration component " + cmpt);

  fits_close_file(ff, &status);
  check_fitsio_status(status);

  return root + "/bcf/" + res;
}


InstPar::InstPar(int tm)
{
  std::string instparfn = lookup_cal("tm"+std::to_string(tm), "GEOM");
  std::printf("Reading INSTPAR file %s\n", instparfn.c_str());

  int status = 0;
  fitsfile* ff;

  fits_open_file(&ff, instparfn.c_str(), READONLY, &status);
  check_fitsio_status(status);

  fits_movnam_hdu(ff, BINARY_TBL, const_cast<char*>("INSTPAR"), 0, &status);
  check_fitsio_status(status);

  // lookup column number and read single value
  auto read_single_col_val = [&ff, &status](const char *colname)
    {
      int col;
      fits_get_colnum(ff, CASEINSEN, const_cast<char*>(colname), &col, &status);
      check_fitsio_status(status);
      double val;
      fits_read_col(ff, TDOUBLE, col, 1, 1, 1, 0, &val, 0, &status);
      check_fitsio_status(status);
      return val;
    };

  x_optax = read_single_col_val("X_OPTAX");
  y_optax = read_single_col_val("Y_OPTAX");
  x_platescale = read_single_col_val("X_PLATESCALE");
  y_platescale = read_single_col_val("Y_PLATESCALE");
  x_ccdpix = read_single_col_val("X_CCDPIX");
  y_ccdpix = read_single_col_val("Y_CCDPIX");
  x_ref = read_single_col_val("X_REF");
  y_ref = read_single_col_val("Y_REF");

  pixscale_x = x_platescale/3600.;
  pixscale_y = y_platescale/3600.;
  inv_pixscale_x = 3600./x_platescale;
  inv_pixscale_y = 3600./y_platescale;

  fits_close_file(ff, &status);
  check_fitsio_status(status);
  std::printf("  - done\n");
}
