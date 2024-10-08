#ifndef PARS_HH
#define PARS_HH

#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "attitude.hh"
#include "deadcor.hh"
#include "detmap.hh"
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
             AttitudeTable,DetMap,
             DeadCorTable> loadEventFile() const;
  InstPar loadInstPar() const;
  Mask loadMask() const;
  std::unique_ptr<ProjMode> createProjMode() const;
  Point imageCentre() const;

  std::vector<std::string> getHeaders() const;

public:
  enum projmodetype {
    AVERAGE_FOV, AVERAGE_FOV_SKY, WHOLE_DET,
    RADIAL, RADIAL_SYM,
    BOX
  };

  enum runmodetype : int { IMAGE, EXPOS };

public:
  // Mode to use
  runmodetype mode;

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

  // bitpix for exposure map
  int bitpix;

  // time delta for exposure map
  double deltat;

  // additionally mask source out to this radius (pixels)
  float masksrcrad;

  // filenames
  std::string evt_fn;
  std::string mask_fn;
  std::string out_fn;
  std::string gti_fn;
};

#endif
