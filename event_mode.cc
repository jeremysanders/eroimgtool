#include <cmath>
#include <cstdio>
#include <mutex>
#include <thread>
#include <vector>

#include <fitsio.h>

#include "event_mode.hh"
#include "common.hh"
#include "geom.hh"
#include "coords.hh"
#include "image.hh"
#include "poly_fill.hh"
#include "events.hh"

// this is similar to image_mode, but we write a fits event table instead

namespace
{
  struct Chunk
  {
    double src_ra, src_dec;
    size_t start, size;
  };

  struct EventOut
  {
    float dx, dy, pi;
  };
}

static void processEvents(std::vector<Chunk>& chunks,
                          std::mutex& mutex,
                          const EventTable& events,
                          Pars pars, GTITable gti, AttitudeTable att,
                          DetMap detmap, Mask mask, InstPar instpar,
                          std::vector<EventOut>& final_out)
{
  auto projmode = pars.createProjMode();
  CoordConv coordconv(instpar);

  // working events
  std::vector<EventOut> evts_out;
  evts_out.reserve(8192);

  for(;;)
    {
      // get next time to process
      Chunk chunk;
      {
        std::lock_guard<std::mutex> lock(mutex);
        if(chunks.empty())
          {
            // merge our events with the total and return
            final_out.insert(final_out.end(), evts_out.begin(), evts_out.end());
            return;
          }
        chunk = chunks.back();
        chunks.pop_back();
      }

      for(size_t i=chunk.start; i!=std::min(chunk.start+chunk.size, events.num_entries); ++i)
        {
          // skip events on bad pixels
          if( detmap.getMap(events.time[i])(events.rawx[i]-1, events.rawy[i]-1) == 0.f )
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
          auto [src_ccdx, src_ccdy] = coordconv.radec2ccd(chunk.src_ra, chunk.src_dec);

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
          relpt = mat.apply(relpt);

          // calculate coordinates in image and add to pixel
          evts_out.push_back({relpt.x, relpt.y, events.pi[i]});
        }

    } // chunks
}


void eventMode(const Pars& pars)
{
  InstPar instpar = pars.loadInstPar();
  auto [events, gti, att, detmap, deadc] = pars.loadEventFile();
  Mask mask = pars.loadMask();

  pars.createProjMode()->message();
  pars.showSources();

  std::printf("Building event list\n");

  // build up a list of sources and chunks of photons
  constexpr size_t chunk_size = 400;
  std::vector<Chunk> chunks;
  for(size_t i=0; i < events.num_entries; i += chunk_size)
    for(auto& srcpos : pars.sources)
      chunks.push_back({srcpos[0], srcpos[1], i, chunk_size});
  // we want it reversed, as vector processing starts from the end
  std::reverse(chunks.begin(), chunks.end());

  std::vector<EventOut> evts_out;
  std::mutex mutex;

  if(pars.threads <= 1)
    {
      processEvents(chunks, mutex,
                    events,
                    pars, gti, att, detmap, mask, instpar, evts_out);
    }
  else
    {
      std::vector<std::thread> threads;
      for(unsigned i=0; i != pars.threads; ++i)
        threads.emplace_back(processEvents,
                             std::ref(chunks), std::ref(mutex),
                             std::ref(events),
                             pars, gti, att, detmap, mask, instpar,
                             std::ref(evts_out));
      for(auto& thread : threads)
        thread.join();
    }


  // make fits file and write data
  std::printf("  - writing output events to %s\n", pars.out_fn.c_str());

  int status = 0;
  fitsfile* ff;
  fits_create_file(&ff, pars.out_fn.c_str(), &status);
  check_fitsio_status(status);

  // make table
  int tfields = 3;
  const char *ttype[] = {"DX", "DY", "PI"};
  const char *tform[] = {"E", "E", "E"};
  const char *tunit[] = {"PIX", "PIX", ""};

  fits_insert_btbl(ff, 0, tfields, const_cast<char**>(ttype), const_cast<char**>(tform),
                   const_cast<char**>(tunit), "EROEVT", 0, &status);
  check_fitsio_status(status);

  // write columns
  std::vector<float> vals;
  for(const auto& e : evts_out)
    vals.push_back(e.dx);
  fits_write_col(ff, TFLOAT, 1, 1, 1, vals.size(), &vals[0], &status);
  vals.clear();
  for(const auto& e : evts_out)
    vals.push_back(e.dy);
  fits_write_col(ff, TFLOAT, 2, 1, 1, vals.size(), &vals[0], &status);
  vals.clear();
  for(const auto& e : evts_out)
    vals.push_back(e.pi);
  fits_write_col(ff, TFLOAT, 3, 1, 1, vals.size(), &vals[0], &status);
  check_fitsio_status(status);

  // close file
  fits_close_file(ff, &status);
  check_fitsio_status(status);
}
