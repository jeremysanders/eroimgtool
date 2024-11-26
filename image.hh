#ifndef IMAGE_HH
#define IMAGE_HH

#include <algorithm>
#include <cassert>
#include <string>
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

  void operator=(T v)
  {
    arr = v;
  }
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

  // Create a new image based on a subset of this image, with width w
  // and height h, from (x0,y0). If the bounds are exceeded, the
  // output image is smaller.
  Image<T> subrect(unsigned x0, unsigned y0, unsigned w, unsigned h)
  {
    const unsigned nw = std::min(w, xw-x0);
    const unsigned nh = std::min(h, yw-y0);
    Image<T> sub(nw, nh);
    for(unsigned y=0; y<nh; ++y)
      for(unsigned x=0; x<nw; ++x)
        sub(x,y) = (*this)(x+x0,y+y0);
    return sub;
  }

  // return image where y and x are swapped
  Image<T> transpose() const
  {
    Image<T> out(yw, xw);
    for(unsigned y=0; y<yw; ++y)
      for(unsigned x=0; x<xw; ++x)
        out(y, x) = (*this)(x, y);
    return out;
  }

  T min() const { return arr.min(); }
  T max() const { return arr.max(); }

  unsigned xw, yw;
  std::valarray<T> arr;
};

void write_fits_image(const std::string& filename,
                      const Image<int>& img,
                      float xc, float yc, float pixscale,
                      bool overwrite=true);

void write_fits_image(const std::string& filename,
                      const Image<float>& img,
                      float xc, float yc, float pixscale,
                      bool overwrite=true,
                      int bitpix=-32);

Image<float> read_fits_image(const std::string& filename);

#endif
