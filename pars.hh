#ifndef PARS_HH
#define PARS_HH

#include <string>
#include <tuple>

#include "attitude.hh"
#include "badpix.hh"
#include "events.hh"
#include "gti.hh"
#include "instpar.hh"
#include "mask.hh"
#include "mode.hh"

class Pars
{
public:
  Pars();
  std::tuple<EventTable,GTITable,
             AttitudeTable,BadPixTable> loadEventFile() const;
  InstPar loadInstPar() const;
  Mask loadMask() const;
  Mode* createMode() const;

public:
  enum modetype {AVERAGE_FOV, AVERAGE_FOV_SKY};

public:
  // TM to process
  int tm;
  // source position
  double src_ra, src_dec;
  // PI range
  float pimin, pimax;

  // operating mode
  modetype mode;

  // output image size
  unsigned xw, yw;
  // output pixel size
  float pixsize;

  // filenames
  std::string evt_fn;
  std::string mask_fn;
  std::string out_fn;
};

#endif
