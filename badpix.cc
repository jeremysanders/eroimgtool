#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include "badpix.hh"
#include "common.hh"
#include "build_poly.hh"

BadPixTable::BadPixTable(fitsfile *ff, int tm)
  : cache_ti(-1)
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

  num_entries = nrows;
  rawx.resize(num_entries);
  rawy.resize(num_entries);
  yextent.resize(num_entries);
  timemin.resize(num_entries);
  timemax.resize(num_entries);

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

  // replace times with inf or -inf as appropriate
  const double inf = std::numeric_limits<double>::infinity();
  for(size_t i=0; i != num_entries; ++i)
    {
      if(!std::isfinite(timemin[i])) timemin[i] = -inf;
      if(!std::isfinite(timemax[i])) timemax[i] = +inf;
    }

  // get list of times where things change
  tedge.push_back(-inf);
  tedge.insert(tedge.end(), timemin.begin(), timemin.end());
  tedge.insert(tedge.end(), timemax.begin(), timemax.end());
  tedge.push_back(+inf);

  std::sort(tedge.begin(), tedge.end());
  tedge.erase( std::unique(tedge.begin(), tedge.end()), tedge.end() );
}

Image<int> BadPixTable::buildMask(double t)
{
  Image<int> mask(384, 384, 1);

  // TODO: account for if neighbouring pixels need to be masked
  for(size_t i=0; i != num_entries; ++i)
    if(t>=timemin[i] && t<timemax[i])
      {
        for(int y=rawy[i]; y<rawy[i]+yextent[i]; ++y)
          mask(rawx[i]-1, y-1) = 0;
      }

  return mask;
}

const PolyVec& BadPixTable::getPolyMask(double t)
{
  if(cache_ti < 0 || t < tedge[cache_ti] || t >= tedge[cache_ti+1])
    {
      // get index for this valid time
      for(size_t i=0; i != tedge.size(); ++i)
        if(t >= tedge[i])
          {
            cache_ti = i;
            break;
          }

      // now build up map
      Image<int> badpixmap(buildMask(t));

      // convert map to polygons
      cache_poly = mask_to_polygons(badpixmap);
    }

  return cache_poly;
}
