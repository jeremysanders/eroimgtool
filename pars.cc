#include <cstdio>
#include <stdexcept>

#include <fitsio.h>

#include "common.hh"
#include "pars.hh"

Pars::Pars()
: tm(1),
  src_ra(0), src_dec(0),
  pimin(300), pimax(2300),
  mode(AVERAGE_FOV),
  xw(512), yw(512),
  pixsize(1)
{
}

std::tuple<EventTable,GTITable,AttitudeTable,BadPixTable>
Pars::loadEventFile() const
{
  int status = 0;
  fitsfile* ff;

  std::printf("Opening event file %s\n", evt_fn.c_str());
  fits_open_file(&ff, evt_fn.c_str(), READONLY, &status);
  check_fitsio_status(status);

  BadPixTable badpix(ff, tm);
  GTITable gti(ff, tm);
  AttitudeTable att(ff, tm);
  EventTable events(ff);
  events.filter_tm(tm);
  events.filter_pi(pimin, pimax);

  fits_close_file(ff, &status);
  check_fitsio_status(status);

  return std::make_tuple(events, gti, att, badpix);
}

InstPar Pars::loadInstPar() const
{
  return InstPar(tm);
}

Mask Pars::loadMask() const
{
  return Mask(mask_fn);
}

Mode* Pars::createMode() const
{
  switch(mode)
    {
    case AVERAGE_FOV:
      return new ModeAverageFoV();
    case AVERAGE_FOV_SKY:
      return new ModeAverageFoVSky();
    default:
      throw std::runtime_error("Invalid mode");
    }
}
