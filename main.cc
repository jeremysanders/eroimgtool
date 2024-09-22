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
#include "mode.hh"

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
  Events events(ff);
  events.filter_tm(tm);
  events.filter_pi(300,2300);
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

  ModeAverageFoV mode;

  for(size_t i=0; i!=events.num_entries; ++i)
    {
      // skip events on bad pixels
      if( bp.getMask(events.time[i])(events.rawx[i]-1, events.rawy[i]-1) == 0 )
        continue;

      Point evtpt(events.ccdx[i], events.ccdy[i]);

      // get attitude at time of event
      auto [att_ra, att_dec, att_roll] = att.interpolate(events.time[i]);
      coordconv.updatePointing(att_ra, att_dec, att_roll);

      // ignore masked regions
      PolyVec ccd_polys(mask.as_ccd_poly(coordconv));
      bool masked = false;
      for(auto& poly : ccd_polys)
        if( poly.is_inside(evtpt) )
          {
            masked = true; break;
          }
      if(masked)
        continue;

      // get ccd coordinates of source
      auto [src_ccdx, src_ccdy] = coordconv.radec2ccd(src_ra, src_dec);

      // skip if source is outsite allowed region
      Point srcccd(src_ccdx, src_ccdy);
      if( ! mode.source_valid(srcccd) )
        continue;

      // compute relative coordinates of photon
      Point origin = mode.origin(srcccd);
      Point relpt = evtpt - origin;

      // apply any necessary rotation for mode
      Point delpt = srcccd - Point(instpar.x_ref, instpar.y_ref);
      auto mat = mode.rotation_matrix(att_roll, delpt);

      relpt = Point(relpt.x*mat[0]+relpt.y*mat[1],
                    relpt.x*mat[2]+relpt.y*mat[3]);

      // calculate coordinates in image and add to pixel
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
