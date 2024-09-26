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
  deltat(0.01)
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

  EventTable events(ff);
  events.filter_tm(tm);
  events.filter_pi(pimin, pimax);

  GTITable gti(ff, tm);
  AttitudeTable att(ff, tm);
  BadPixTable badpix(ff, tm);
  DeadCorTable deadc(ff, tm);

  fits_close_file(ff, &status);
  check_fitsio_status(status);

  return std::make_tuple(events, gti, att, badpix, deadc);
}

InstPar Pars::loadInstPar() const
{
  return InstPar(tm);
}

Mask Pars::loadMask() const
{
  return Mask(mask_fn);
}

ProjMode* Pars::createProjMode() const
{
  switch(projmode)
    {
    case AVERAGE_FOV:
      return new ProjModeAverageFoV();
    case AVERAGE_FOV_SKY:
      return new ProjModeAverageFoVSky();
    case WHOLE_DET:
      return new ProjModeDet();
    default:
      throw std::runtime_error("Invalid mode");
    }
}

Point Pars::imageCentre() const
{
  return Point(xw/2, yw/2);
}
