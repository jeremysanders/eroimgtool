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

  fits_movnam_hdu(ff, BINARY_TBL, const_cast<char*>("CIF"), 0, &status);
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

  int c_x_optax, c_y_optax, c_x_platescale, c_y_platescale,
    c_x_ccdpix, c_y_ccdpix, c_x_ref, c_y_ref, c_timelag;
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("X_OPTAX"),
                  &c_x_optax, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("Y_OPTAX"),
                  &c_y_optax, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("X_PLATESCALE"),
                  &c_x_platescale, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("Y_PLATESCALE"),
                  &c_y_platescale, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("X_CCDPIX"),
                  &c_x_ccdpix, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("Y_CCDPIX"),
                  &c_y_ccdpix, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("X_REF"),
                  &c_x_ref, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("Y_REF"),
                  &c_y_ref, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("TIMELAG"),
                  &c_timelag, &status);
  check_fitsio_status(status);

  fits_read_col(ff, TDOUBLE, c_x_optax, 1, 1, 1, 0, &x_optax, 0, &status);
  fits_read_col(ff, TDOUBLE, c_y_optax, 1, 1, 1, 0, &y_optax, 0, &status);
  fits_read_col(ff, TDOUBLE, c_x_platescale, 1, 1, 1, 0, &x_platescale, 0, &status);
  fits_read_col(ff, TDOUBLE, c_y_platescale, 1, 1, 1, 0, &y_platescale, 0, &status);
  fits_read_col(ff, TDOUBLE, c_x_ccdpix, 1, 1, 1, 0, &x_ccdpix, 0, &status);
  fits_read_col(ff, TDOUBLE, c_y_ccdpix, 1, 1, 1, 0, &y_ccdpix, 0, &status);
  fits_read_col(ff, TDOUBLE, c_x_ref, 1, 1, 1, 0, &x_ref, 0, &status);
  fits_read_col(ff, TDOUBLE, c_y_ref, 1, 1, 1, 0, &y_ref, 0, &status);
  fits_read_col(ff, TDOUBLE, c_timelag, 1, 1, 1, 0, &timelag, 0, &status);
  check_fitsio_status(status);

  fits_close_file(ff, &status);
  check_fitsio_status(status);
  std::printf("  Done\n");
}
