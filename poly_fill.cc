#include <algorithm>
#include <cmath>
#include <vector>

#include "poly_fill.hh"
#include "common.hh"

// swap a and b, if a>b
// this uses min and max to make it branchless
inline void order_pair(float& a, float& b)
{
  float t=a;
  a = std::min(a,b);
  b = std::max(t,b);
}

static void sortSmall(std::vector<float>& vals)
{
  // these are sort networks, optimized for small numbers
  // https://bertdobbelaere.github.io/sorting_networks.html
  switch(vals.size())
    {
    case 0:
    case 1:
      break;
    case 2:
      order_pair(vals[0], vals[1]);
      break;
    case 3:
      order_pair(vals[0], vals[2]);
      order_pair(vals[0], vals[1]);
      order_pair(vals[1], vals[2]);
      break;
    case 4:
      order_pair(vals[0], vals[2]);
      order_pair(vals[1], vals[3]);
      order_pair(vals[0], vals[1]);
      order_pair(vals[2], vals[3]);
      order_pair(vals[1], vals[2]);
      break;
    case 5:
      order_pair(vals[0], vals[3]);
      order_pair(vals[1], vals[4]);
      order_pair(vals[0], vals[2]);
      order_pair(vals[1], vals[3]);
      order_pair(vals[0], vals[1]);
      order_pair(vals[2], vals[4]);
      order_pair(vals[1], vals[2]);
      order_pair(vals[3], vals[4]);
      order_pair(vals[2], vals[3]);
      break;
    case 6:
      order_pair(vals[0], vals[5]);
      order_pair(vals[1], vals[3]);
      order_pair(vals[2], vals[4]);
      order_pair(vals[1], vals[2]);
      order_pair(vals[3], vals[4]);
      order_pair(vals[0], vals[3]);
      order_pair(vals[2], vals[5]);
      order_pair(vals[0], vals[1]);
      order_pair(vals[2], vals[3]);
      order_pair(vals[4], vals[5]);
      order_pair(vals[1], vals[2]);
      order_pair(vals[3], vals[4]);
      break;
    default:
      std::sort(vals.begin(), vals.end());
    }
}

// alterantive to doing (i+1)%npts without division
inline int nextwrap(int i, int npts)
{
  return (i+1==npts) ? 0 : i+1;
}

void fillPoly(const Poly& poly, Image<uint8_t>& outimg, uint8_t val)
{
  const int xw = outimg.xw;
  const int yw = outimg.yw;
  const int npts = poly.size();

  const Rect bounds = poly.bounds();
  const int ylo = std::max(int(std::floor(bounds.tl.y)), 0);
  const int yhi = std::min(int(std::ceil(bounds.br.y)), yw-1);

  // compute gradients (hopefully won't use the infinite ones, as the
  // if statement below shouldn't be used)
  std::vector<float> grads(npts);
  for(int i=0; i<npts; ++i)
    grads[i] = (poly[nextwrap(i,npts)].x - poly[i].x) /
      (poly[nextwrap(i,npts)].y - poly[i].y);

  std::vector<float> xs;
  xs.reserve(8);

  for(int y=ylo; y<=yhi; ++y)
    {
      const float yf = y;

      // find which polygons cross this y
      for(int pi=0; pi<npts; ++pi)
        {
          const float x1 = poly[pi].x;
          const float y1 = poly[pi].y;
          const float y2 = poly[nextwrap(pi,npts)].y;

          // have we crossed the line?
          if( ((y1<=yf) && (y2>yf)) || ((y1>yf) && (y2<=yf)) )
            {
              // solve for x where line crosses
              xs.emplace_back(x1 + grads[pi]*(yf-y1));
            }
        }
      sortSmall(xs);

      // fill in pixels between pairs of points where it cross a line
      const int nxs = xs.size();
      for(int i=0; i+1 < nxs; i+= 2)
        {
          const int xlo = std::max(int(std::ceil(xs[i])), 0);
          const int xhi = std::min(int(std::floor(xs[i+1])), xw-1);
          // little optimization as we're accessing pixels on a line
          for(int x=xlo; x<=xhi; ++x)
            outimg.arr[y*xw+x] = val;
        }
      xs.clear();
    }
}

/*
int main()
{
  int xw=256;
  int yw=256;
  Image<uint8_t> img(xw, yw);

  Poly poly1;
  poly1.add(Point(-50,50));
  poly1.add(Point(50,100));
  poly1.add(Point(100,100));
  poly1.add(Point(120,100));
  poly1.add(Point(120,120));
  poly1.add(Point(100,120));
  poly1.add(Point(100,100));
  poly1.add(Point(50,100));
  poly1.add(Point(50,200));
  poly1.add(Point(200,200));
  poly1.add(Point(200,50));

  poly1.rotate(15*DEG2RAD);
  fillPoly(poly1, img, 1);

  Image<int> imcpy(xw, yw);
  for(int y=0; y<yw; ++y)
    for(int x=0; x<xw; ++x)
      imcpy(x, y) = img(x, y);

  write_fits_image("out.fits", imcpy, xw/2, yw/2, 1);

  return 0;
}

*/
