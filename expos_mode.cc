#include <cmath>
#include <cstdio>
#include <memory>
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
  TimeSeg() {}
  TimeSeg(size_t _idx, double _t, float _dt) : idx(_idx), t(_t), dt(_dt) {}

  size_t idx;
  double t;
  double dt;
};

static void processGTIs(size_t num,
                        std::vector<TimeSeg>& times,
                        std::mutex& mutex,
                        Pars pars, GTITable gti, AttitudeTable att,
                        BadPixTable bp,
                        Mask mask, InstPar instpar,
                        Image<double>& finalimg)
{
  std::unique_ptr<ProjMode> projmode(pars.createProjMode());
  CoordConv coordconv(instpar);
  Point imgcen = pars.imageCentre();

  // output image
  Image<double> img(pars.xw, pars.yw, 0.f);
  // masked image during time step
  Image<uint8_t> imgt(pars.xw, pars.yw, 0);

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
      auto [src_ccdx, src_ccdy] = coordconv.radec2ccd(pars.src_ra, pars.src_dec);

      // skip if source is outsite allowed region
      Point srcccd(src_ccdx, src_ccdy);
      if( ! projmode->sourceValid(srcccd) )
        continue;

      if( timeseg.idx % 200 == 0 )
        std::printf("Iteration %.1f%% (t=%.1f)\n", timeseg.idx*100./num, timeseg.t);

      Point delpt = srcccd - Point(instpar.x_ref, instpar.y_ref);
      auto mat = projmode->rotationMatrix(att_roll, delpt);
      Point projorigin = projmode->origin(srcccd);

      // polygons defining detector
      PolyVec detpolys(bp.getPolyMask(timeseg.t));
      applyShiftRotationScaleShift(detpolys, mat, projorigin, 1/pars.pixsize, imgcen);

      // polygons with bad regions
      PolyVec maskedpolys(mask.as_ccd_poly(coordconv));
      applyShiftRotationScaleShift(maskedpolys, mat, projorigin, 1/pars.pixsize, imgcen);

      imgt = 0;
      for(auto& poly: detpolys)
        fillPoly(poly, imgt, 1);
      for(auto& poly: maskedpolys)
        fillPoly(poly, imgt, 0);

      int npix = img.xw * img.yw;
      for(int i=0; i<npix; ++i)
        img.arr[i] += (imgt.arr[i]==1) ? timeseg.dt : 0.;

      //fillPoly2(detpolys, maskedpolys, img, timeseg.dt);

    } // input times

}

void exposMode(const Pars& pars)
{
  InstPar instpar = pars.loadInstPar();
  auto [events, gti, att, bp, deadc] = pars.loadEventFile();

  Mask mask = pars.loadMask();
  //mask.writeRegion("test.reg");

  std::unique_ptr<ProjMode> projmode(pars.createProjMode());
  projmode->message();

  std::printf("Building exposure map\n");

  // put times in vector
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
          timesegs.emplace_back(timesegs.size(), t, deltat*deadcf);
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
                  pars, gti, att, bp,
                  mask, instpar, sumimg);
    }
  else
    {
      std::vector<std::thread> threads;
      for(unsigned i=0; i != pars.threads; ++i)
        threads.emplace_back(processGTIs,
                             num, std::ref(timesegs), std::ref(mutex),
                             pars, gti, att, bp, mask, instpar,
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
  write_fits_image(pars.out_fn, writeimg, imgcen.x, imgcen.y, pars.pixsize);
}
