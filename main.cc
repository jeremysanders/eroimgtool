#include <cmath>
#include <cstdio>

#include <fitsio.h>

#include "badpix.hh"
#include "gti.hh"
#include "attitude.hh"
#include "events.hh"
#include "common.hh"
#include "geom.hh"
#include "coords.hh"
#include "image.hh"
#include "instpar.hh"
#include "mask.hh"

int main()
{
  const int tm = 2;

  InstPar instpar(tm);

  int status = 0;

  fitsfile* ff;

  const char *filename = "em01_056102_020_ML00001_004_c946/evt.fits.gz";
  std::printf("Opening file %s\n", filename);
  fits_open_file(&ff, filename, READONLY, &status);
  check_fitsio_status(status);

  BadPixTable bp(ff, tm);
  GTITable gti(ff, tm);
  AttitudeTable att(ff, tm);
  Events events(ff, tm);
  Mask mask("em01_056102_020_ML00001_004_c946/030_mask_final.fits.gz");

  std::printf("Making image\n");
  double src_ra = 57.3466206;
  double src_dec = -11.9909090;

  CoordConv coordconv(instpar.x_platescale, instpar.y_platescale,
                      instpar.x_ref, instpar.y_ref);

  Image<int> outimg(512, 512);
  float xc = 255;
  float yc = 255;
  float pixscale = 1; //1/8.f;

  for(size_t i=0; i!=events.num_entries; ++i)
    {
      // skip events on bad pixels
      if( bp.getMask(events.time[i])(events.rawx[i]-1, events.rawy[i]-1) == 0 )
        continue;

      auto [att_ra, att_dec, att_roll] = att.interpolate(events.time[i]);
      coordconv.updatePointing(att_ra, att_dec, att_roll);

      // ignore masked regions
      PolyVec ccd_polys(mask.as_ccd_poly(coordconv));
      bool ok = true;
      Point pt(events.ccdx[i], events.ccdy[i]);
      for(auto& poly : ccd_polys)
        {
          ok = ok && !poly.is_inside(pt);
          if(!ok) break;
        }

      if(!ok)
        continue;

      auto [src_ccdx, src_ccdy] = coordconv.radec2ccd(src_ra, src_dec);
      // src_ccdx = 192.5; src_ccdy = 192.5;

      Point relpt = pt - Point(src_ccdx, src_ccdy);
      Point scalept = relpt/pixscale + Point(xc, yc);
      int px = int(std::round(scalept.x));
      int py = int(std::round(scalept.y));
      if(px>=0 && px<int(outimg.xw) && py>=0 && py<int(outimg.yw))
        outimg(px, py) += 1;
    }

  write_fits_image("test.fits", outimg, xc, yc, pixscale);


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
