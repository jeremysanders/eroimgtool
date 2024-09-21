#ifndef COMMON_HH
#define COMMON_HH

#include <algorithm>
#include <numeric>

#include <fitsio.h>

// throw an exception if status!=0 with fitsio error
void check_fitsio_status(int status);

// quick helper to read a whole fits column
void read_fits_column(fitsfile *ff, const char* name,
                      int type, long nrows, void *retn);

// move to hdu with name given
void move_fits_hdu(fitsfile* ff, const char* name);

// maths
constexpr double PI = 3.14159265358979323846264338327950288;
constexpr double DEG2RAD = PI / 180.;
constexpr double RAD2DEG = 180. / PI;

// erosita constants
constexpr unsigned CCD_XW = 384;
constexpr unsigned CCD_YW = 384;

template<class T> T sqr(T a)
{
  return a*a;
}

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

// produce a list of indices which sort
template <class T> std::vector<size_t> argsort(const std::vector<T>& v)
{
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

// return entries in v, indexed by idx
// i.e. v[idx[0]], v[idx[1]], ...
template <class T> std::vector<T> selidx(const std::vector<T>& v,
                                         const std::vector<size_t>& idx)
{
  std::vector<T> retn;
  retn.reserve(idx.size());
  for(size_t i : idx)
    retn.push_back(v[i]);
  return retn;
}

#endif
