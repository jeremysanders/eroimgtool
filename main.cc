#include <cmath>
#include <cstdio>
#include <memory>

#include "common.hh"
#include "geom.hh"
#include "coords.hh"
#include "image.hh"
#include "pars.hh"

void imageMode(const Pars& pars)
{
  InstPar instpar = pars.loadInstPar();
  auto [events, gti, att, bp] = pars.loadEventFile();
  Mask mask = pars.loadMask();

  std::printf("Building image\n");

  CoordConv coordconv(instpar);

  Image<int> outimg(pars.xw, pars.yw);
  Point ptc = pars.imageCentre();

  std::unique_ptr<ProjMode> projmode(pars.createProjMode());

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
      PolyVec ccd_maskedpolys(mask.as_ccd_poly(coordconv));
      if( is_inside(ccd_maskedpolys, evtpt) )
        continue;

      // get ccd coordinates of source
      auto [src_ccdx, src_ccdy] = coordconv.radec2ccd(pars.src_ra, pars.src_dec);

      // skip if source is outsite allowed region
      Point srcccd(src_ccdx, src_ccdy);
      if( ! projmode->sourceValid(srcccd) )
        continue;

      // compute relative coordinates of photon
      Point origin = projmode->origin(srcccd);
      Point relpt = evtpt - origin;

      // apply any necessary rotation for mode
      Point delpt = srcccd - Point(instpar.x_ref, instpar.y_ref);
      auto mat = projmode->rotationMatrix(att_roll, delpt);
      relpt = mat.rotate(relpt);

      // calculate coordinates in image and add to pixel
      Point scalept = relpt/pars.pixsize + ptc;
      int px = int(std::round(scalept.x));
      int py = int(std::round(scalept.y));
      if(px>=0 && px<int(outimg.xw) && py>=0 && py<int(outimg.yw))
        outimg(px, py) += 1;
    }

  std::printf("  - writing output image to %s\n", pars.out_fn.c_str());
  write_fits_image(pars.out_fn, outimg, ptc.x, ptc.y, pars.pixsize);
}

void exposMode(const Pars& pars)
{
  InstPar instpar = pars.loadInstPar();
  auto [events, gti, att, bp] = pars.loadEventFile();
  Mask mask = pars.loadMask();

  std::printf("Building exposure map\n");

  CoordConv coordconv(instpar);

  Image<float> outimg(pars.xw, pars.yw, 0.f);
  Point ptc = pars.imageCentre();

  std::unique_ptr<ProjMode> mode(pars.createProjMode());

  std::printf("  - writing output image to %s\n", pars.out_fn.c_str());
  write_fits_image(pars.out_fn, outimg, ptc.x, ptc.y, pars.pixsize);

/*
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
*/
}

int main()
{
  Pars pars;
  pars.tm = 2;
  pars.evt_fn = "em01_056102_020_ML00001_004_c946/evt.fits.gz";
  pars.mask_fn = "em01_056102_020_ML00001_004_c946/030_mask_final.fits.gz";
  pars.out_fn = "test.fits";
  pars.src_ra = 57.3466206;
  pars.src_dec = -11.9909090;
  pars.pixsize = 1;

  //imageMode(pars);
  exposMode(pars);

  return 0;
}
