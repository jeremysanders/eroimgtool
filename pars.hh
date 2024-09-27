#ifndef PARS_HH
#define PARS_HH

#include <string>
#include <tuple>
#include <vector>

#include "attitude.hh"
#include "badpix.hh"
#include "deadcor.hh"
#include "events.hh"
#include "geom.hh"
#include "gti.hh"
#include "instpar.hh"
#include "mask.hh"
#include "proj_mode.hh"

class Pars
{
public:
  Pars();
  std::tuple<EventTable,GTITable,
             AttitudeTable,BadPixTable,
             DeadCorTable> loadEventFile() const;
  InstPar loadInstPar() const;
  Mask loadMask() const;
  ProjMode* createProjMode() const;
  Point imageCentre() const;

public:
  enum projmodetype {AVERAGE_FOV, AVERAGE_FOV_SKY, WHOLE_DET, RADIAL};

public:
  // TM to process
  int tm;
  // source position
  double src_ra, src_dec;
  // PI range
  float pimin, pimax;

  // operating mode
  projmodetype projmode;

  // arguments to projection mode
  std::vector<float> projargs;

  // number of threads to use
  unsigned threads;

  // output image size
  unsigned xw, yw;
  // output pixel size
  float pixsize;

  // time delta for exposure map
  double deltat;

  // additionally mask source out to this radius (pixels)
  float masksrcrad;

  // filenames
  std::string evt_fn;
  std::string mask_fn;
  std::string out_fn;
};

#endif
