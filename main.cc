#include <cmath>
#include <cstdio>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>

#include "common.hh"
#include "geom.hh"
#include "coords.hh"
#include "image.hh"
#include "pars.hh"
#include "poly_fill.hh"

#include "image_mode.hh"
#include "expos_mode.hh"

int main()
{
  Pars pars;
  pars.tm = 2;
  pars.evt_fn = "em01_056102_020_ML00001_004_c946/evt.fits.gz";
  pars.mask_fn = "em01_056102_020_ML00001_004_c946/030_mask_final.fits.gz";
  pars.out_fn = "test.fits";
  pars.src_ra = 57.3466206;
  pars.src_dec = -11.9909090;
  pars.pixsize = 1/5.f;///8.f;
  pars.threads = 16;

  pars.xw = 1024;
  pars.yw = 1024;

  // pars.xw *= 4;
  // pars.yw *= 4;
  // pars.pixsize /= 4;

  //pars.deltat = 0.1;
  //pars.projmode = Pars::WHOLE_DET;
  //pars.projmode = Pars::AVERAGE_FOV_SKY;

  imageMode(pars);
  //exposMode(pars);

  return 0;
}
