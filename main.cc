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

  BadPixTable bp(ff, 2);
  GTITable gti(ff, 2);
  AttitudeTable att(ff, 2);

  double t = 6.41406e+08;
  PolyVec polys = bp.getPolyMask(t);

  Image<float> img(512,512,0.f);

  std::printf("a\n");
  for(auto& bppoly : polys)
    {
      Rect bound = bppoly.bounds();

      for(int y=0; y<int(img.yw); ++y)
        for(int x=0; x<int(img.xw); ++x)
          {
            float xf = x/512.*400.;
            float yf = y/512.*400.;

            if( (xf+1<bound.tl.x) || (yf+1<bound.tl.y) || (xf>bound.br.x) ||
                (yf>bound.br.y) )
              continue;

            Poly pixp;
            pixp.add(Point(xf,yf));
            pixp.add(Point(xf,yf+1));
            pixp.add(Point(xf+1,yf+1));
            pixp.add(Point(xf+1,yf));

            Poly clipped = poly_clip(bppoly, pixp);
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
