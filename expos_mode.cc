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

// min and max of 4 values
static inline float min4(float a, float b, float c, float d)
{
  return std::min(std::min(a, b), std::min(c, d));
}
static inline float max4(float a, float b, float c, float d)
{
  return std::max(std::max(a, b), std::max(c, d));
}

// do rectangles defined from x1->x2 and y1->y2 (incl) overlap?
static inline bool rectoverlap(int ax1, int ax2, int ay1, int ay2,
                               int bx1, int bx2, int by1, int by2)
{
  return ax1<=bx2 && ax2>=bx1 && ay2>=by1 && ay1<=by2;
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
  Image<double> img(pars.xw, pars.yw, 0.);

  // image during time step
  Image<float> imgt(pars.xw, pars.yw, 0.f);

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

      Point srcccd(src_ccdx, src_ccdy);
      Point delpt = srcccd - Point(instpar.x_ref, instpar.y_ref);
      Point projorigin = projmode->origin(srcccd);

      // matrix to go from detector -> image, including pixel size
      auto mat = projmode->rotationMatrix(att_roll, delpt);
      mat.scale(1/pars.pixsize);

      // matrix to go from image -> detector, including pixel size
      auto matrev = projmode->rotationMatrix(-att_roll, delpt);
      matrev.scale(pars.pixsize);

      // find coordinates of detector corners in output image
      const Point ic1 = mat.apply(Point(0,0) - projorigin) + imgcen;
      const Point ic2 = mat.apply(Point(CCD_XW,0) - projorigin) + imgcen;
      const Point ic3 = mat.apply(Point(0,CCD_YW) - projorigin) + imgcen;
      const Point ic4 = mat.apply(Point(CCD_XW,CCD_YW) - projorigin) + imgcen;

      // calculate integer bounding box of detector in output image
      const int ic_xlo = int(std::floor(min4(ic1.x, ic2.x, ic3.x, ic4.x)));
      const int ic_xhi = int(std::ceil (max4(ic1.x, ic2.x, ic3.x, ic4.x)));
      const int ic_ylo = int(std::floor(min4(ic1.y, ic2.y, ic3.y, ic4.y)));
      const int ic_yhi = int(std::ceil (max4(ic1.y, ic2.y, ic3.y, ic4.y)));

      // skip if there's no overlap between detector and output image
      if(! rectoverlap(ic_xlo, ic_xhi, ic_ylo, ic_yhi, -1, pars.xw, -1, pars.yw))
        continue;

      if( timeseg.idx % 500 == 0 )
        std::printf("Iteration %5.1f%% (t=%.1f)\n", timeseg.idx*100./num, timeseg.t);

      // detector map for time
      const Image<float>& dmimg = detmap.getMap(timeseg.t);

      // these are the ranges to iterate over
      const int minx = std::clamp(ic_xlo-1, 0, int(pars.xw)-1);
      const int maxx = std::clamp(ic_xhi+1, 0, int(pars.xw)-1);
      const int miny = std::clamp(ic_ylo-1, 0, int(pars.yw)-1);
      const int maxy = std::clamp(ic_yhi+1, 0, int(pars.yw)-1);

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

static std::vector<TimeSeg> applySampling(const std::vector<TimeSeg>& timesegs, int samples)
{
  std::printf("  - making %d samples in time\n", samples);

  // get total time in segments
  double tott = 0;
  for(size_t i=0; i != timesegs.size(); ++i)
    tott += double(timesegs[i].dt);
  const double deltat = tott / samples; // size of new segments

  // evenly sample along cumulative time
  size_t ts = 0; // current time segment
  double tsum = double(timesegs[ts].dt); // total time up to current time segment

  std::vector<TimeSeg> newsegs;
  for(int i=0; i<samples; ++i)
    {
      // cumulative time of new segment
      const double t = (i+0.5) * deltat;

      // find corresponding old segment in cumulative time
      while(t>tsum)
        {
          ++ts;
          tsum += double(timesegs[ts].dt);
        }

      // create new segment
      newsegs.emplace_back( TimeSeg({
            timesegs[ts].src_ra,
            timesegs[ts].src_dec,
            size_t(i),
            timesegs[ts].t,
            float(deltat)
            }) );
    }
  return newsegs;
}

void exposMode(const Pars& pars)
{
  InstPar instpar = pars.loadInstPar();
  auto [events, gti, att, detmap, deadc] = pars.loadEventFile();

  Mask mask = pars.loadMask();

  auto projmode = pars.createProjMode();
  projmode->message();
  pars.showSources();

  std::printf("Building exposure map\n");

  CoordConv coordconv(instpar);

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

          auto [att_ra, att_dec, att_roll] = att.interpolate(t);
          coordconv.updatePointing(att_ra, att_dec, att_roll);

          float deadcf = deadc.interpolate(t);
          for(auto& srcpos : pars.sources)
            {
              auto [src_ccdx, src_ccdy] = coordconv.radec2ccd(srcpos[0], srcpos[1]);

              // add time if source is inside region
              Point srcccd(src_ccdx, src_ccdy);
              if( projmode->sourceValid(srcccd) )
                {
                  timesegs.emplace_back( TimeSeg({
                        srcpos[0], srcpos[1],
                        timesegs.size(), t, float(deltat*deadcf)}) );
                }
            } // sources
        } // times in GTI
    } // GTIs

  // sample if requested, but only if it results in fewer calculations
  if(pars.samples > 0 && pars.samples < int(timesegs.size()))
    {
      timesegs = applySampling(timesegs, pars.samples);
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
