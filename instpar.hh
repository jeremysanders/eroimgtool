#ifndef INSTPAR_HH
#define INSTPAR_HH

#include <string>

std::string lookup_cal(const std::string& subdir, const std::string& cmpt);

class InstPar
{
public:
  InstPar(int tm);

  double x_optax, y_optax;
  double x_platescale, y_platescale;
  double x_ccdpix, y_ccdpix;
  double x_ref, y_ref;

  double pixscale_x, pixscale_y;
  double inv_pixscale_x, inv_pixscale_y;
};


#endif
