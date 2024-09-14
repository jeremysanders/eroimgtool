#ifndef COORDS_HH
#define COORDS_HH

#include <cmath>
#include <tuple>

constexpr double DEG2RAD = M_PI / 180.;
constexpr double RAD2DEG = 180. / M_PI;

class CoordConv
{
public:

  CoordConv(double _x_platescale, double _y_platescale,
            double _x_ref, double _y_ref)
    : x_platescale(_x_platescale), y_platescale(_y_platescale),
      x_ref(_x_ref), y_ref(_y_ref),
      rad2xpix(1/(x_platescale * (DEG2RAD / 3600.))),
      rad2ypix(1/(y_platescale * (DEG2RAD / 3600.))),
      ra0(0), sindec0(0), cosdec0(0), rsin(0), rcos(0)
    {
    }

  // set telescope pointing from attitude
  void updatePointing(double _ra0, double _dec0, double _roll0)
  {
    ra0 = _ra0;
    sindec0 =  std::sin(_dec0*DEG2RAD);
    cosdec0 =  std::cos(_dec0*DEG2RAD);
    double rtheta = (_roll0-90.)*DEG2RAD;
    rsin = std::sin(rtheta);
    rcos = std::cos(rtheta);
  }

  // convert RA, Dec to CCD coordinates
  std::tuple<double, double> radec2ccd(double ra, double dec) const
  {
    double diffra = (ra - ra0)*DEG2RAD;
    double dsinra = std::sin(diffra);
    double dcosra = std::cos(diffra);

    double sindec = std::sin(dec*DEG2RAD);
    double cosdec = std::cos(dec*DEG2RAD);

    double d1s =    dsinra * cosdec;
    double dh  =  - cosdec * dcosra;
    double d1c =    sindec0 * sindec - dh * cosdec0;
    double dx  =    std::atan2(d1s,d1c);

    double d2s = dh * sindec0 + cosdec0 * sindec;
    double d2c = std::sqrt(1.0 - d2s * d2s);

    double dy  = -std::atan2(d2s,d2c);

    // now rotate about roll
    double rx = dx*rcos - dy*rsin;
    double ry = dx*rsin + dy*rcos;

    // apply plate scale
    double ccdx = rx*rad2xpix + x_ref;
    double ccdy = ry*rad2ypix + y_ref;

    return std::tuple<double,double>(ccdx, ccdy);
  }

private:
  
  double x_platescale, y_platescale, x_ref, y_ref, rad2xpix, rad2ypix;
  double ra0, sindec0, cosdec0, rsin, rcos;
};


#endif
