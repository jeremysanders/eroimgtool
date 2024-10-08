#include <cstdio>
#include <stdexcept>

#include <fitsio.h>

#include "common.hh"
#include "pars.hh"

Pars::Pars()
: mode(IMAGE),
  tm(1),
  src_ra(0), src_dec(0),
  pimin(300), pimax(2300),
  projmode(AVERAGE_FOV),
  threads(1),
  xw(512), yw(512),
  pixsize(1),
  bitpix(-32),
  deltat(0.01),
  masksrcrad(0)
{
}

std::tuple<EventTable,GTITable,AttitudeTable,BadPixTable,DeadCorTable>
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
  BadPixTable badpix(ff, tm);
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

  return std::make_tuple(events, gti, att, badpix, deadc);
}

InstPar Pars::loadInstPar() const
{
  return InstPar(tm);
}

Mask Pars::loadMask() const
{
  Mask mask(mask_fn);
  if(masksrcrad>0)
    mask.setSrcMask(src_ra, src_dec, masksrcrad);
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

std::vector<std::string> Pars::getHeaders() const
{
  std::vector<std::string> hdrs;

  hdrs.emplace_back("--tm=" + std::to_string(tm));
  hdrs.emplace_back("--ra=" + std::to_string(src_ra));
  hdrs.emplace_back("--dec=" + std::to_string(src_dec));
  hdrs.emplace_back("--pi-min=" + std::to_string(pimin));
  hdrs.emplace_back("--pi-max=" + std::to_string(pimax));
  hdrs.emplace_back("--proj=" + std::to_string(projmode));

  if(!projargs.empty())
    {
      hdrs.emplace_back("--proj-args");
      std::string temp;
      for(auto v : projargs)
        temp = temp + std::to_string(v) + " ";
      hdrs.emplace_back(temp);
    }

  hdrs.emplace_back("--threads=" + std::to_string(threads));
  hdrs.emplace_back("--xw=" + std::to_string(xw));
  hdrs.emplace_back("--yw=" + std::to_string(yw));
  hdrs.emplace_back("--pixsize=" + std::to_string(pixsize));
  hdrs.emplace_back("--delta-t=" + std::to_string(deltat));
  hdrs.emplace_back("--mask-src-rad=" + std::to_string(masksrcrad));

  if(!mask_fn.empty())
    hdrs.emplace_back("--mask=" + mask_fn);
  if(!gti_fn.empty())
    hdrs.emplace_back("--gti=" + gti_fn);

  hdrs.emplace_back(std::to_string(mode));
  hdrs.emplace_back(evt_fn);
  hdrs.emplace_back(out_fn);

  return hdrs;
}
