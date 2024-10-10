#include <cmath>
#include <cstdio>
#include <mutex>
#include <stdexcept>
#include <thread>

#include "expos_mode.hh"
#include "common.hh"
#include "geom.hh"
#include "coords.hh"
#include "image.hh"
#include "poly_fill.hh"

struct TimeSeg
{
  double src_ra, src_dec;
  size_t idx;
  double t;
  float dt;
};

// find minimum, subtract one, then limit to range minval to maxval
inline int min4intclamp(float a, float b, float c, float d, int minval, int maxval)
{
  int x = int( std::floor(std::min(std::min(a, b), std::min(c, d))) );
  return std::clamp(x-1, minval, maxval);
}

// find maximum, add one, then limit to range minval to maxval
inline int max4intclamp(float a, float b, float c, float d, int minval, int maxval)
{
  int x = int( std::ceil(std::max(std::max(a, b), std::max(c, d))) );
  return std::clamp(x+1, minval, maxval);
}

static void processGTIs(size_t num,
                        std::vector<TimeSeg>& times,
                        std::mutex& mutex,
                        Pars pars, GTITable gti, AttitudeTable att,
                        DetMap detmap,
                        Mask mask, InstPar instpar,
                        Image<double>& finalimg)
{
  auto projmode = pars.createProjMode();
  CoordConv coordconv(instpar);
  Point imgcen = pars.imageCentre();

  // output image
  Image<double> img(pars.xw, pars.yw, 0.f);

  // image during time step
  Image<float> imgt(pars.xw, pars.yw, 0);

  for(;;)
    {
      // get next time to process
      TimeSeg timeseg;
      {
        std::lock_guard<std::mutex> lock(mutex);
        if(times.empty())
          {
            // add our part to the total and return
            finalimg.arr += img.arr;
            return;
          }
        timeseg = times.back();
        times.pop_back();
      }
      auto [att_ra, att_dec, att_roll] = att.interpolate(timeseg.t);
      coordconv.updatePointing(att_ra, att_dec, att_roll);

      // get ccd coordinates of source
      auto [src_ccdx, src_ccdy] = coordconv.radec2ccd(timeseg.src_ra, timeseg.src_dec);

      // skip if source is outsite allowed region
      Point srcccd(src_ccdx, src_ccdy);
      if( ! projmode->sourceValid(srcccd) )
        continue;

      if( timeseg.idx % 500 == 0 )
        std::printf("Iteration %5.1f%% (t=%.1f)\n", timeseg.idx*100./num, timeseg.t);

      Point delpt = srcccd - Point(instpar.x_ref, instpar.y_ref);
      Point projorigin = projmode->origin(srcccd);

      // detector map for time
      const Image<float>& dmimg = detmap.getMap(timeseg.t);

      // matrix to go from detector -> image, including pixel size
      auto mat = projmode->rotationMatrix(att_roll, delpt);
      mat.scale(1/pars.pixsize);

      // matrix to go from image -> detector, including pixel size
      auto matrev = projmode->rotationMatrix(-att_roll, delpt);
      matrev.scale(pars.pixsize);

      // find range of coordinates of detector in output image
      Point ic1 = mat.apply(Point(0,0) - projorigin) + imgcen;
      Point ic2 = mat.apply(Point(CCD_XW,0) - projorigin) + imgcen;
      Point ic3 = mat.apply(Point(0,CCD_YW) - projorigin) + imgcen;
      Point ic4 = mat.apply(Point(CCD_XW,CCD_YW) - projorigin) + imgcen;

      const int minx = min4intclamp(ic1.x, ic2.x, ic3.x, ic4.x, 0, pars.xw-1);
      const int maxx = max4intclamp(ic1.x, ic2.x, ic3.x, ic4.x, 0, pars.xw-1);
      const int miny = min4intclamp(ic1.y, ic2.y, ic3.y, ic4.y, 0, pars.yw-1);
      const int maxy = max4intclamp(ic1.y, ic2.y, ic3.y, ic4.y, 0, pars.yw-1);

      // iterate over output image, pixel by pixel
      float* optr = &imgt.arr[0];
      for(int i=0; i<miny*int(pars.xw); ++i)
        *optr++ = 0.f;
      for(int y=miny; y<=maxy; ++y)
        {
          for(int x=0; x<minx; ++x)
            *optr++ = 0.f;
          for(int x=minx; x<=maxx; ++x)
            {
              // rotate around imgcen and move to origin
              Point det = matrev.apply(Point(x,y)-imgcen) + projorigin;

              // the funny +16. -16 is to ensure correct rounding around zero
              // without this -0.5 is rounded 0 due to truncation
              // could use int(std::floor()) instead, but is quite a lot slower
              int dix = int(det.x+(16.f-0.5f))-16;
              int diy = int(det.y+(16.f-0.5f))-16;

              if(dix>=0 && diy>=0 && dix<int(CCD_XW) && diy<int(CCD_YW))
                *optr++ = dmimg(dix, diy);
              else
                *optr++ = 0.f;
            }
          for(int x=maxx+1; x<int(pars.xw); ++x)
            *optr++ = 0.f;
        }
      for(int i=(maxy+1)*int(pars.xw); i<int(pars.xw*pars.yw); ++i)
        *optr++ = 0.f;

      // zero out polygons with bad regions
      PolyVec maskedpolys(mask.as_ccd_poly(coordconv));
      applyShiftRotationShift(maskedpolys, mat, projorigin, imgcen);
      for(auto& poly: maskedpolys)
        fillPoly(poly, imgt, 0);

      int npix = img.xw * img.yw;
      for(int i=0; i<npix; ++i)
        img.arr[i] += imgt.arr[i] * timeseg.dt;

    } // input times

}

void exposMode(const Pars& pars)
{
  InstPar instpar = pars.loadInstPar();
  auto [events, gti, att, detmap, deadc] = pars.loadEventFile();

  Mask mask = pars.loadMask();

  pars.createProjMode()->message();
  pars.showSources();

  std::printf("Building exposure map\n");

  // put sources and times in vector
  std::vector<TimeSeg> timesegs;
  for(int gtii=0; gtii<int(gti.num); ++gtii)
    {
      double tstart = gti.start[gtii];
      double tstop = gti.stop[gtii];

      if(tstop<=tstart)
        throw std::runtime_error("invalid GTI found");
      int numt = int(std::ceil((tstop - tstart) / pars.deltat));
      double deltat = (tstop - tstart) / numt;

      for(int ti=0; ti<numt; ++ti)
        {
          double t = tstart + (ti+0.5)*deltat;
          float deadcf = deadc.interpolate(t);
          for(auto& srcpos : pars.sources)
            timesegs.emplace_back( TimeSeg({
                  srcpos[0], srcpos[1], timesegs.size(), t, float(deltat*deadcf)}) );
        }
    }
  // we process them in order of time, so reverse so we pop from back
  std::reverse(timesegs.begin(), timesegs.end());

  //timesegs.clear();
  //timesegs.emplace_back(0, 635459963.8, 0.1);

  // summed output image
  Image<double> sumimg(pars.xw, pars.yw, 0.f);

  std::mutex mutex;

  size_t num = timesegs.size();
  if(pars.threads <= 1)
    {
      processGTIs(num, timesegs, mutex,
                  pars, gti, att, detmap,
                  mask, instpar, sumimg);
    }
  else
    {
      std::vector<std::thread> threads;
      for(unsigned i=0; i != pars.threads; ++i)
        threads.emplace_back(processGTIs,
                             num, std::ref(timesegs), std::ref(mutex),
                             pars, gti, att, detmap, mask, instpar,
                             std::ref(sumimg));
      for(auto& thread : threads)
        thread.join();
    }

  Image<float> writeimg(pars.xw, pars.yw);
  for(unsigned y=0; y<pars.yw; ++y)
    for(unsigned x=0; x<pars.xw; ++x)
      writeimg(x,y) = float(sumimg(x,y));

  Point imgcen = pars.imageCentre();
  std::printf("  - writing output image to %s\n", pars.out_fn.c_str());
  write_fits_image(pars.out_fn, writeimg, imgcen.x, imgcen.y, pars.pixsize, true, pars.bitpix);
}
