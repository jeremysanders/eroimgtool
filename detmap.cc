#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

#include "detmap.hh"
#include "common.hh"

DetMap::DetMap(fitsfile *ff, int tm)
  : cache_ti(-1),
    cache_map(CCD_XW, CCD_YW)
{
  int status = 0;

  std::string hduname = std::string("BADPIX") + std::to_string(tm);
  std::printf("  - Opening bad pixel extension %s\n", hduname.c_str());
  move_fits_hdu(ff, hduname.c_str());

  long nrows;
  fits_get_num_rows(ff, &nrows, &status);
  check_fitsio_status(status);

  num_entries = nrows;
  rawx.resize(num_entries);
  rawy.resize(num_entries);
  yextent.resize(num_entries);
  timemin.resize(num_entries);
  timemax.resize(num_entries);

  read_fits_column(ff, "RAWX", TINT, nrows, &rawx[0]);
  read_fits_column(ff, "RAWY", TINT, nrows, &rawy[0]);
  read_fits_column(ff, "YEXTENT", TINT, nrows, &yextent[0]);
  read_fits_column(ff, "TIMEMIN", TDOUBLE, nrows, &timemin[0]);
  read_fits_column(ff, "TIMEMAX", TDOUBLE, nrows, &timemax[0]);

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

void DetMap::checkCache(double t)
{
  if(cache_ti < 0 || t < tedge[cache_ti] || t >= tedge[cache_ti+1])
    {
      size_t i=0;
      while( i+1 < tedge.size() && tedge[i+1]<t )
        ++i;
      cache_ti = i;

      buildMapImage(t);
    }
}

void DetMap::buildMapImage(double t)
{
  cache_map = 1.f;

  for(size_t i=0; i != num_entries; ++i)
    if(t>=timemin[i] && t<timemax[i])
      {
        // expand bad pixel by 1
        int ylo = std::max(rawy[i]-1, 1);
        int yhi = std::min(rawy[i]+yextent[i]+1-1, int(CCD_YW));
        int xlo = std::max(rawx[i]-1, 1);
        int xhi = std::min(rawx[i]+1, int(CCD_XW));

        for(int y=ylo; y<=yhi; ++y)
          for(int x=xlo; x<=xhi; ++x)
            cache_map(x-1, y-1) = 0.f;
      }

  for(int i=0; i<384; i++)
    {
      cache_map(0, i) = 0.f;
      cache_map(383, i) = 0.f;
      cache_map(i, 0) = 0.f;
      cache_map(i, 383) = 0.f;
    }
}
