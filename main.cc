#include <fitsio.h>
#include <cstdio>
#include "badpix.hh"
#include "gti.hh"
#include "attitude.hh"
#include "common.hh"
#include "geom.hh"
#include "coords.hh"

#include "image.hh"

int main()
{
  Poly poly;
  poly.add(Point(1,0));
  poly.add(Point(0,1));
  poly.rotate(90*DEG2RAD);
  for(auto p: poly.pts)
    std::printf("%g %g\n", p.x, p.y);

  int status = 0;

  fitsfile* ff;

  const char *filename = "em01_108141_020_ML00003_003_c946/evt.fits.gz";
  std::printf("Opening file %s\n", filename);
  fits_open_file(&ff, filename, READONLY, &status);
  check_fitsio_status(status);

  BadPixTable bp(ff, 2);
  GTITable gti(ff, 2);
  AttitudeTable att(ff, 2);

  fits_close_file(ff, &status);
  check_fitsio_status(status);

  return 0;
}
