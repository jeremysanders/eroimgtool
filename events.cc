#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <limits>

#include "events.hh"
#include "common.hh"

EventTable::EventTable(fitsfile *ff)
{
  int status = 0;

  std::printf("  - Opening EVENTS extension\n");
  move_fits_hdu(ff, "EVENTS");

  long nrows;
  fits_get_num_rows(ff, &nrows, &status);

  check_fitsio_status(status);

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

  // combine rawx/y and subx/y
  ccdx.reserve(nrows); ccdy.reserve(nrows);
  for(long i=0; i<nrows; ++i)
    {
      ccdx.push_back(rawx[i]+subx[i]);
      ccdy.push_back(rawy[i]+suby[i]);
    }

  // sort entries by time (normally sorted anyway)
  std::vector<size_t> sort_idx = argsort(time);
  do_filter(sort_idx);

  std::printf("    - successfully read %ld entries\n", nrows);
}

void EventTable::filter_tm(int tm)
{
  std::vector<size_t> idxs;
  for(size_t i=0; i != rawx.size(); ++i)
    if(tm_nr[i] == tm)
      idxs.push_back(i);

  do_filter(idxs);
  std::printf("    - filtered to TM%d, giving %ld entries\n", tm, num_entries);
}

void EventTable::filter_pi(float pimin, float pimax)
{
  std::vector<size_t> idxs;
  for(size_t i=0; i != rawx.size(); ++i)
    if(pi[i]>=pimin && pi[i]<pimax)
      idxs.push_back(i);

  do_filter(idxs);
  std::printf("    - filtered from PI=%g:%g, giving %ld entries\n",
              pimin, pimax, num_entries);
}

// filter all columns to have indices given
void EventTable::do_filter(const std::vector<size_t>& sel)
{
  rawx = selidx(rawx, sel);
  rawy = selidx(rawy, sel);
  tm_nr = selidx(tm_nr, sel);
  ra = selidx(ra, sel);
  dec = selidx(dec, sel);
  time = selidx(time, sel);
  pi = selidx(pi, sel);
  subx = selidx(subx, sel);
  suby = selidx(suby, sel);
  ccdx = selidx(ccdx, sel);
  ccdy = selidx(ccdy, sel);
  num_entries = rawx.size();
}
