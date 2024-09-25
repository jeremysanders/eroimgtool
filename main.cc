#include <cmath>
#include <cstdio>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>

#include "common.hh"
#include "geom.hh"
#include "coords.hh"
#include "image.hh"
#include "pars.hh"
#include "poly_fill.hh"

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

static std::vector<Rect> getPolysBounds(const PolyVec& polys)
{
  std::vector<Rect> bounds;
  bounds.reserve(polys.size());
  for(auto& poly : polys)
    bounds.emplace_back( poly.bounds() );
  return bounds;
}

struct TimeSeg
{
  TimeSeg() {}
  TimeSeg(size_t _idx, double _t, float _dt) : idx(_idx), t(_t), dt(_dt) {}

  size_t idx;
  double t;
  double dt;
};

void processGTIs(size_t num,
                 std::vector<TimeSeg>& times,
                 std::mutex& mutex,
                 Pars pars, GTITable gti, AttitudeTable att, BadPixTable bp,
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
        //img.arr[i] += timeseg.dt * imgt.arr[i];
        img.arr[i] += (imgt.arr[i]==1) ? timeseg.dt : 0.;
        //if(imgt.arr[i]) img.arr[i] += timeseg.dt;
        //img.arr[i] += imgt.arr[i] ? timeseg.dt : 0.;

      //fillPoly2(detpolys, maskedpolys, img, timeseg.dt);

    } // input times

}

// 1:16 - tern
// 1:21 - mult
// 1:27 - if
// 1:20 - tern2

void exposMode(const Pars& pars)
{
  InstPar instpar = pars.loadInstPar();
  auto [events, gti, att, bp] = pars.loadEventFile();

  Mask mask = pars.loadMask();
  //mask.simplify();
  //mask.writeRegion("test.reg");

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
          timesegs.emplace_back(timesegs.size(), t, deltat);
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
                  pars, gti, att, bp, mask, instpar, sumimg);
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

int main()
{
  Pars pars;
  pars.tm = 2;
  pars.evt_fn = "em01_056102_020_ML00001_004_c946/evt.fits.gz";
  pars.mask_fn = "em01_056102_020_ML00001_004_c946/030_mask_final.fits.gz";
  pars.out_fn = "test.fits";
  pars.src_ra = 57.3466206;
  pars.src_dec = -11.9909090;
  pars.pixsize = 1/5.f;///8.f;
  pars.threads = 4;

  pars.xw = 1024;
  pars.yw = 1024;

  // pars.xw *= 4;
  // pars.yw *= 4;
  // pars.pixsize /= 4;

  //pars.deltat = 0.1;
  //pars.projmode = Pars::WHOLE_DET;
  //pars.projmode = Pars::AVERAGE_FOV_SKY;

  //imageMode(pars);
  exposMode(pars);

  return 0;
}
