#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>

#include "events.hh"
#include "common.hh"

// filter vector according to tm
template <class T> static void filter_vector_tm
(std::vector<T>& vals, const std::vector<short>& tm_nr, short tm)
{
  assert(vals.size() == tm_nr.size());

  size_t j = 0;
  for(size_t i=0; i != tm_nr.size(); ++i)
    {
      if(tm_nr[i] == tm)
        {
          vals[j] = vals[i];
          ++j;
        }
    }

  vals.resize(j);
}

Events::Events(fitsfile *ff, int tm)
{
  int status = 0;

  std::printf("  Opening EVENTS extension\n");
  fits_movnam_hdu(ff, BINARY_TBL, const_cast<char*>("EVENTS"),
                  0, &status);
  check_fitsio_status(status);

  // get indices of columns
  int c_rawx, c_rawy, c_tm_nr, c_ra, c_dec, c_time,
    c_pi, c_subx, c_suby;
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("RAWX"),
                  &c_rawx, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("RAWY"),
                  &c_rawy, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("TM_NR"),
                  &c_tm_nr, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("RA"),
                  &c_ra, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("DEC"),
                  &c_dec, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("TIME"),
                  &c_time, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("PI"),
                  &c_pi, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("SUBX"),
                  &c_subx, &status);
  fits_get_colnum(ff, CASEINSEN, const_cast<char*>("SUBY"),
                  &c_suby, &status);

  long nrows;
  fits_get_num_rows(ff, &nrows, &status);

  check_fitsio_status(status);

  std::vector<short> tm_nr;
  rawx.resize(nrows);
  rawy.resize(nrows);
  tm_nr.resize(nrows);
  ra.resize(nrows);
  dec.resize(nrows);
  time.resize(nrows);
  pi.resize(nrows);
  subx.resize(nrows);
  suby.resize(nrows);

  // read in table
  fits_read_col(ff, TSHORT, c_rawx, 1, 1, nrows, 0, &rawx[0], 0, &status);
  fits_read_col(ff, TSHORT, c_rawy, 1, 1, nrows, 0, &rawy[0], 0, &status);
  fits_read_col(ff, TSHORT, c_tm_nr, 1, 1, nrows, 0, &tm_nr[0], 0, &status);
  fits_read_col(ff, TDOUBLE, c_ra, 1, 1, nrows, 0, &ra[0], 0, &status);
  fits_read_col(ff, TDOUBLE, c_dec, 1, 1, nrows, 0, &dec[0], 0, &status);
  fits_read_col(ff, TDOUBLE, c_time, 1, 1, nrows, 0, &time[0], 0, &status);
  fits_read_col(ff, TFLOAT, c_pi, 1, 1, nrows, 0, &pi[0], 0, &status);
  fits_read_col(ff, TFLOAT, c_subx, 1, 1, nrows, 0, &subx[0], 0, &status);
  fits_read_col(ff, TFLOAT, c_suby, 1, 1, nrows, 0, &suby[0], 0, &status);

  check_fitsio_status(status);

  std::printf("    - successfully read %ld entries\n", nrows);

  // filter rows to just get events for TM
  filter_vector_tm(rawx, tm_nr, tm);
  filter_vector_tm(rawy, tm_nr, tm);
  filter_vector_tm(ra, tm_nr, tm);
  filter_vector_tm(dec, tm_nr, tm);
  filter_vector_tm(time, tm_nr, tm);
  filter_vector_tm(pi, tm_nr, tm);
  filter_vector_tm(subx, tm_nr, tm);
  filter_vector_tm(suby, tm_nr, tm);

  // sort entries by time (normally sorted anyway)
  std::vector<size_t> sort_idx = argsort(time);
  rawx = selidx(rawx, sort_idx);
  rawy = selidx(rawy, sort_idx);
  ra = selidx(ra, sort_idx);
  dec = selidx(dec, sort_idx);
  time = selidx(time, sort_idx);
  pi = selidx(pi, sort_idx);
  subx = selidx(subx, sort_idx);
  suby = selidx(suby, sort_idx);

  num_entries = rawx.size();

  // combine rawx/y and subx/y
  for(size_t i=0; i!=num_entries; ++i)
    {
      ccdx.push_back(rawx[i]+subx[i]);
      ccdy.push_back(rawy[i]+suby[i]);
    }

  std::printf("    - filtered to TM%d to %ld entries\n", tm, num_entries);
}
