#include <cstdio>
#include <stdexcept>

#include <fitsio.h>

#include "common.hh"
#include "pars.hh"

Pars::Pars()
: tm(1),
  src_ra(0), src_dec(0),
  pimin(300), pimax(2300),
  projmode(AVERAGE_FOV),
  threads(1),
  xw(512), yw(512),
  pixsize(1),
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
