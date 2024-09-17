#ifndef COMMON_HH
#define COMMON_HH

// throw an exception if status!=0 with fitsio error
void check_fitsio_status(int status);

// maths
constexpr double PI = 3.14159265358979323846264338327950288;
constexpr double DEG2RAD = PI / 180.;
constexpr double RAD2DEG = 180. / PI;

// erosita constants
constexpr unsigned CCD_XW = 384;
constexpr unsigned CCD_YW = 384;

// divide, rounding up
template<class T> constexpr T div_round_up(T a, T b)
{
  return a/b + (a%b != 0);
}

// clip to range
template <typename T> T clip(T v, T minv, T maxv)
{
  return v < minv ? minv : (v > maxv ? maxv : v);
}

#endif
