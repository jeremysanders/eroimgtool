#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include "badpix.hh"
#include "common.hh"
#include "build_poly.hh"

BadPixTable::BadPixTable(fitsfile *ff, int tm)
  : cache_ti(-1),
    cache_mask(CCD_XW, CCD_YW)
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

void BadPixTable::buildMask(double t)
{
  cache_mask = 1;

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
            cache_mask(x-1, y-1) = 0;
      }
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

      // we subdivide the image into tiles, as we can more easily
      // check bounds for a small polygon and the polygon is much
      // simpler if there are bad pixels in the middle
      cache_poly.clear();

      buildMask(t);

      constexpr unsigned tile_size = 32;
      constexpr unsigned num_tile_x = div_round_up(CCD_XW, tile_size);
      constexpr unsigned num_tile_y = div_round_up(CCD_YW, tile_size);
      for(unsigned yt=0;yt<num_tile_y;++yt)
        for(unsigned xt=0;xt<num_tile_x;++xt)
          {
            // chop input
            auto rect = cache_mask.subrect(xt*tile_size, yt*tile_size,
                                           tile_size, tile_size);
            // make polygons
            auto polys = mask_to_polygons(rect);

            // add to cached list of polygons
            const Point p0(xt*tile_size+0.5f, yt*tile_size+0.5f);
            for(auto& poly : polys)
              {
                // if(xt==3 && yt==799)
                //   {
                //     for(int y=0;y<tile_size;y++) {
                //       for(int x=0;x<tile_size;++x) {
                //         std::printf("%c", rect(x,y)>0 ? '#' : '.');
                //       }
                //       std::printf("\n");
                //     }
                //     for(auto pt : poly.pts) {
                //       std::printf("%g %g\n", pt.x, pt.y);
                //     }
                //   }
                // std::printf("%d, %d, %ld\n", xt, yt, poly.size());
                cache_poly.emplace_back(poly + p0);
              }
          }

    }

  return cache_poly;
}
