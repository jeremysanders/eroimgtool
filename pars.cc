#include <cstdio>
#include <stdexcept>

#include <fitsio.h>

#include "common.hh"
#include "pars.hh"

Pars::Pars()
: mode(IMAGE),
  tm(1),
  pimin(300), pimax(2300),
  projmode(AVERAGE_FOV),
  detmapmask(false),
  shadowmask(false),
  threads(1),
  xw(512), yw(512),
  pixsize(1),
  bitpix(-32),
  deltat(0.01),
  samples(-1)
{
}

std::tuple<EventTable,GTITable,AttitudeTable,DetMap,DeadCorTable>
Pars::loadEventFile() const
{
  int status = 0;
  fitsfile* ff;

  std::printf("Opening event file %s\n", evt_fn.c_str());
  fits_open_file(&ff, evt_fn.c_str(), READONLY, &status);
  check_fitsio_status(status);

  GTITable gti(ff, tm);

  EventTable events(ff);
  events.filter_tm(tm);
  events.filter_pi(pimin, pimax);
  events.filter_gti(gti);

  AttitudeTable att(ff, tm);
  DetMap detmap(tm, detmapmask, shadowmask);
  detmap.read(ff);
  DeadCorTable deadc(ff, tm);

  fits_close_file(ff, &status);
  check_fitsio_status(status);

  if(!gti_fn.empty())
    {
      std::printf("Opening GTI file %s\n", gti_fn.c_str());
      fits_open_file(&ff, gti_fn.c_str(), READONLY, &status);
      check_fitsio_status(status);
      GTITable gti2(ff, tm);
      fits_close_file(ff, &status);
      check_fitsio_status(status);

      gti &= gti2;
      std::printf("  - merged GTIs to make %ld elements\n", gti.num);
      events.filter_gti(gti);
    }

  if(!bpix_fn.empty())
    {
      detmap.read(bpix_fn);
    }

  return std::make_tuple(events, gti, att, detmap, deadc);
}

void Pars::showSources() const
{
  std::printf("Defined %ld source(s)\n", sources.size());
  for(auto& src : sources)
    std::printf("  - source (%g,%g)\n", src[0], src[1]);
}

InstPar Pars::loadInstPar() const
{
  return InstPar(tm);
}

Mask Pars::loadMask() const
{
  Mask mask(mask_fn);
  mask.setMaskPts(maskpts);
  return mask;
}

std::unique_ptr<ProjMode> Pars::createProjMode() const
{
  switch(projmode)
    {
    case AVERAGE_FOV:
      return std::make_unique<ProjModeAverageFoV>();
    case AVERAGE_FOV_SKY:
      return std::make_unique<ProjModeAverageFoVSky>();
    case AVERAGE_FULL:
      return std::make_unique<ProjModeAverageFull>();
    case WHOLE_DET:
      return std::make_unique<ProjModeDet>();
    case RADIAL:
      return std::make_unique<ProjModeRadial>(projargs);
    case RADIAL_SYM:
      return std::make_unique<ProjModeRadialSym>(projargs);
    case BOX:
      return std::make_unique<ProjModeBox>(projargs);
    default:
      throw std::runtime_error("Invalid mode");
    }
}

Point Pars::imageCentre() const
{
  return Point(xw/2, yw/2);
}

namespace
{
  template<typename T> std::string str_list(const std::vector<T>& vals)
  {
    std::string retn;
    for(auto& v : vals)
      {
        if(!retn.empty())
          retn += ',';
        retn += std::to_string(v);
      }
    return retn;
  }

  template<typename T, size_t N> std::string str_list(const std::vector<std::array<T,N>>& vals)
  {
    std::string retn;
    for(auto& arr : vals)
      {
        if(! retn.empty())
          retn += ' ';
        std::string sub;
        for(auto& v : arr)
          {
            if(!sub.empty())
              sub += ',';
            sub += std::to_string(v);
          }
        retn += sub;
      }
    return retn;
  }

}

std::vector<std::string> Pars::getHeaders() const
{
  std::vector<std::string> hdrs;

  hdrs.emplace_back("--tm=" + std::to_string(tm));
  if(!sources.empty())
    hdrs.emplace_back("--sources " + str_list(sources));

  hdrs.emplace_back("--pi-min=" + std::to_string(pimin));
  hdrs.emplace_back("--pi-max=" + std::to_string(pimax));
  hdrs.emplace_back("--proj=" + std::to_string(projmode));

  if(!projargs.empty())
    hdrs.emplace_back("--proj-args " + str_list(projargs));

  hdrs.emplace_back("--threads=" + std::to_string(threads));
  hdrs.emplace_back("--xw=" + std::to_string(xw));
  hdrs.emplace_back("--yw=" + std::to_string(yw));
  hdrs.emplace_back("--pixsize=" + std::to_string(pixsize));
  hdrs.emplace_back("--delta-t=" + std::to_string(deltat));

  if(!mask_fn.empty())
    hdrs.emplace_back("--mask=" + mask_fn);
  if(!maskpts.empty())
    hdrs.emplace_back("--mask-pts " + str_list(maskpts));

  if(!gti_fn.empty())
    hdrs.emplace_back("--gti=" + gti_fn);

  hdrs.emplace_back(std::to_string(mode));
  hdrs.emplace_back(evt_fn);
  hdrs.emplace_back(out_fn);

  return hdrs;
}
