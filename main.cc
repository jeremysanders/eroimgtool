#include <fitsio.h>
#include <cstdio>
#include "badpix.hh"
#include "gti.hh"
#include "attitude.hh"
#include "common.hh"
#include "geom.hh"
#include "coords.hh"

#include "image.hh"

int main()
{
  int status = 0;

  fitsfile* ff;

  const char *filename = "em01_108141_020_ML00003_003_c946/evt.fits.gz";
  std::printf("Opening file %s\n", filename);
  fits_open_file(&ff, filename, READONLY, &status);
  check_fitsio_status(status);

  BadPixTable bp(ff, 4);
  GTITable gti(ff, 2);
  AttitudeTable att(ff, 2);

  double t = 6.41406e+08;
  PolyVec polys = bp.getPolyMask(t);

  Image<float> img(512,512,0.f);

  std::printf("a\n");
  for(auto& poly : polys)
    {
      Poly bppoly = poly *1.5;
      bppoly.rotate(20*DEG2RAD);

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
            pixp[0].x = x;   pixp[0].y = y;
            pixp[1].x = x;   pixp[1].y = y+1;
            pixp[2].x = x+1; pixp[2].y = y+1;
            pixp[3].x = x+1; pixp[3].y = y;

            poly_clip(bppoly, pixp, clipped);
            img(x,y) += clipped.area();
          }
    }
  std::printf("b\n");

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

  fits_close_file(ff, &status);
  check_fitsio_status(status);

  return 0;
}
