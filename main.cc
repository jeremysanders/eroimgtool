#include <fitsio.h>
#include <cstdio>
#include "badpix.hh"
#include "gti.hh"
#include "attitude.hh"
#include "events.hh"
#include "common.hh"
#include "geom.hh"
#include "coords.hh"

#include "image.hh"
#include "instpar.hh"

int main()
{
  int status = 0;

  fitsfile* ff;

  const char *filename = "em01_056102_020_ML00001_004_c946/evt.fits.gz";
  std::printf("Opening file %s\n", filename);
  fits_open_file(&ff, filename, READONLY, &status);
  check_fitsio_status(status);

  BadPixTable bp(ff, 2);
  GTITable gti(ff, 2);
  AttitudeTable att(ff, 2);
  Events events(ff, 2);
  InstPar instpar(2);
  return 0;

  double src_ra = 57.3469733;
  double src_dec = -11.9910211;

  for(size_t i=0; i!=events.num_entries; ++i)
    {
      const auto [att_ra, att_dec, att_roll] = att.interpolate(events.time[i]);




    }



  /*
  Image<float> img(384,384,0.f);

  for(auto& poly : polys)
    {
      Poly bppoly = poly;
      Rect bound = bppoly.bounds();

      const int ylo = std::max(0, int(std::floor(bound.tl.y)));
      const int xlo = std::max(0, int(std::floor(bound.tl.x)));
      const int yhi = std::min(int(img.yw)-1, int(std::ceil(bound.br.y)));
      const int xhi = std::min(int(img.xw)-1, int(std::ceil(bound.br.x)));

      // predefine polygon for pixel
      Poly pixp, clipped;
      for(int i=0; i<4; ++i)
        pixp.add(Point());

      for(int y=ylo; y<=yhi; ++y)
        for(int x=xlo; x<=xhi; ++x)
          {
            pixp[0].x = x+0.5f; pixp[0].y = y+0.5f;
            pixp[1].x = x+0.5f; pixp[1].y = y+1.5f;
            pixp[2].x = x+1.5f; pixp[2].y = y+1.5f;
            pixp[3].x = x+1.5f; pixp[3].y = y+0.5f;

            poly_clip(bppoly, pixp, clipped);
            img(x,y) += clipped.area();
          }
    }

  FILE* f = std::fopen("test.dat", "w");
  for(int y=0; y<int(img.yw); ++y)
    {
      for(int x=0; x<int(img.xw); ++x)
        {
          std::fprintf(f, "%.4f ", img(x,y));
        }
      std::fprintf(f, "\n");
    }
  std::fclose(f);
  */

  fits_close_file(ff, &status);
  check_fitsio_status(status);

  return 0;
}
