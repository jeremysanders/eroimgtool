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

  std::printf("  - Opening EVENTS extension\n");
  move_fits_hdu(ff, "EVENTS");

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
  read_fits_column(ff, "RAWX", TSHORT, nrows, &rawx[0]);
  read_fits_column(ff, "RAWY", TSHORT, nrows, &rawy[0]);
  read_fits_column(ff, "TM_NR", TSHORT, nrows, &tm_nr[0]);
  read_fits_column(ff, "RA", TDOUBLE, nrows, &ra[0]);
  read_fits_column(ff, "DEC", TDOUBLE, nrows, &dec[0]);
  read_fits_column(ff, "TIME", TDOUBLE, nrows, &time[0]);
  read_fits_column(ff, "PI", TFLOAT, nrows, &pi[0]);
  read_fits_column(ff, "SUBX", TFLOAT, nrows, &subx[0]);
  read_fits_column(ff, "SUBY", TFLOAT, nrows, &suby[0]);

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
