#include <cmath>
#include <cstdio>
#include <mutex>
#include <thread>

#include "image_mode.hh"
#include "common.hh"
#include "geom.hh"
#include "coords.hh"
#include "image.hh"
#include "poly_fill.hh"
#include "events.hh"

static void processEvents(size_t chunk_size,
                          std::vector<size_t>& chunks,
                          std::mutex& mutex,
                          const EventTable& events,
                          Pars pars, GTITable gti, AttitudeTable att, BadPixTable bp,
                          Mask mask, InstPar instpar,
                          Image<int>& finalimg)
{
  auto projmode = pars.createProjMode();
  CoordConv coordconv(instpar);
  Point imgcen = pars.imageCentre();

  // working image
  Image<int> img(pars.xw, pars.yw, 0);

  for(;;)
    {
      // get next time to process
      unsigned evtidx;
      {
        std::lock_guard<std::mutex> lock(mutex);
        if(chunks.empty())
          {
            // add our part to the total and return
            finalimg.arr += img.arr;
            return;
          }
        evtidx = chunks.back();
        chunks.pop_back();
      }

      for(size_t i=evtidx; i!=std::min(evtidx+chunk_size, events.num_entries); ++i)
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
          Point scalept = relpt/pars.pixsize + imgcen;
          int px = int(std::round(scalept.x));
          int py = int(std::round(scalept.y));
          if(px>=0 && px<int(img.xw) && py>=0 && py<int(img.yw))
            img(px, py) += 1;
        }

    } // chunks
}


void imageMode(const Pars& pars)
{
  InstPar instpar = pars.loadInstPar();
  auto [events, gti, att, bp, deadc] = pars.loadEventFile();
  Mask mask = pars.loadMask();

  pars.createProjMode()->message();

  std::printf("Building image\n");

  constexpr size_t chunk_size = 400;
  std::vector<size_t> chunks;
  for(size_t i=0; i < events.num_entries; i += chunk_size)
    chunks.push_back(i);
  std::reverse(chunks.begin(), chunks.end());

  Image<int> sumimg(pars.xw, pars.yw, 0);

  std::mutex mutex;

  if(pars.threads <= 1)
    {
      processEvents(chunk_size, chunks, mutex,
                    events,
                    pars, gti, att, bp, mask, instpar, sumimg);
    }
  else
    {
      std::vector<std::thread> threads;
      for(unsigned i=0; i != pars.threads; ++i)
        threads.emplace_back(processEvents,
                             chunk_size, std::ref(chunks), std::ref(mutex),
                             std::ref(events),
                             pars, gti, att, bp, mask, instpar,
                             std::ref(sumimg));
      for(auto& thread : threads)
        thread.join();
    }


  std::printf("  - writing output image to %s\n", pars.out_fn.c_str());
  Point imgcen = pars.imageCentre();
  write_fits_image(pars.out_fn, sumimg, imgcen.x, imgcen.y, pars.pixsize);
}
