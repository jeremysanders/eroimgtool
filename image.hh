#ifndef IMAGE_HH
#define IMAGE_HH

#include <cassert>
#include <valarray>
#include <vector>

#include "geom.hh"

// thin 2D wrapper to array
template <class T> class Image
{
public:
  Image(unsigned _xw, unsigned _yw)
    : xw(_xw), yw(_yw), arr(xw*yw) {}
  Image(unsigned _xw, unsigned _yw, T val)
    : xw(_xw), yw(_yw), arr(val, xw*yw) {}

  T operator()(unsigned x, unsigned y) const
  {
    assert(x<xw && y<yw);
    return arr[x+y*xw];
  }
  T& operator()(unsigned x, unsigned y)
  {
    assert(x<xw && y<yw);
    return arr[x+y*xw];
  }
  unsigned size() const { return xw*yw; }

  unsigned xw, yw;
  std::valarray<T> arr;
};

std::vector<Poly> mask_poly(const Image<int>& img);

#endif
