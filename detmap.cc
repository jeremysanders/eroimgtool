#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

#include "detmap.hh"
#include "common.hh"
#include "instpar.hh"

DetMap::DetMap(fitsfile *ff, int tm, bool detmapmask, bool shadowmask)
  : cache_ti(-1),
    init_map(CCD_XW, CCD_YW),
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

  // setup standard detector map, etc
  if(detmapmask)
    {
      readDetmapMask(tm);
    }
  else
    {
      init_map = 1.f;
    }

  // remove edges initially
  for(unsigned y=0; y<CCD_YW; ++y)
    {
      init_map(0, y) = 0.f;
      init_map(CCD_XW-1, y) = 0.f;
    }
  for(unsigned x=0; x<CCD_XW; ++x)
    {
      init_map(x, 0) = 0.f;
      init_map(x, CCD_YW-1) = 0.f;
    }

  // remove shadowed area from readout if requested
  if(shadowmask)
    {
      for(unsigned y=0; y<15; ++y)
        for(unsigned x=0; x<CCD_XW; ++x)
          init_map(x,y) = 0.f;
    }
}

void DetMap::readDetmapMask(int tm)
{
  std::string fn = lookup_cal("tm"+std::to_string(tm), "DETMAP");
  std::printf("  - Opening DETMAP file %s\n", fn.c_str());

  Image<float> map = read_fits_image(fn);
  if(map.xw != CCD_XW || map.yw != CCD_YW)
    throw std::runtime_error("Invalid detector map size");

  init_map = map;
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
  cache_map = init_map;

  for(size_t i=0; i != num_entries; ++i)
    if(t>=timemin[i] && t<timemax[i])
      {
        int ylo = rawy[i]-1;
        int yhi = rawy[i]+yextent[i]-1;
        int x = rawx[i]-1;

        for(int y=ylo; y<=yhi; ++y)
          {
            // zero pixel
            cache_map(x, y) = 0.f;
            // zero surrounding pixels +-1 in x or y (not both)
            if(x-1>=0)
              cache_map(x-1,y) = 0.f;
            if(x+1<int(CCD_XW))
              cache_map(x+1,y) = 0.f;
            if(y-1>=0)
              cache_map(x,y-1) = 0.f;
            if(y+1<int(CCD_YW))
              cache_map(x,y+1) = 0.f;
          }
      }
}
