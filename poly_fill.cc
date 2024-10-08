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

void fillPoly(const Poly& poly, Image<float>& outimg, float val)
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

namespace
{

  // get gradients between each line segment
  std::vector<float> getGrads(const Poly& poly)
  {
    int npts = poly.size();
    std::vector<float> grads(npts);
    for(int i=0; i<npts; ++i)
      {
        grads[i] = (poly[nextwrap(i,npts)].x - poly[i].x) /
          (poly[nextwrap(i,npts)].y - poly[i].y);
      }
    return grads;
  }

  struct intersect
  {
    intersect() {};
    intersect(float _x, bool _det) : x(_x), det(_det) {};
    float x;
    bool det;
  };

  inline void addCrossings(std::vector<intersect>& isects, float yf, const Poly& poly,
                           const std::vector<float> grads, bool det)
  {
    int npts = poly.size();
    for(int i=0; i<npts; ++i)
      {
        const float x1 = poly[i].x;
        const float y1 = poly[i].y;
        const float y2 = poly[nextwrap(i,npts)].y;

        // have we crossed the line? - if so solve and store
        if( ((y1<=yf) && (y2>yf)) || ((y1>yf) && (y2<=yf)) )
          isects.emplace_back(x1 + grads[i]*(yf-y1), det);
      }
  }

  inline void doFill(int y, Image<float>& outimg, float x1, float x2, float fillval)
  {
    float cx1 = std::ceil(x1);
    float fx2 = std::floor(x2);
    int ix1 = int(cx1);
    int ix2 = int(fx2);

    // draw main part
    for(int x=ix1; x<=ix2; ++x)
      outimg.arr[y*outimg.xw + x] += fillval;
  }
}

void fillPoly2(const PolyVec& detpoly, const PolyVec& maskpoly,
               Image<float>& outimg, float fillval)
{
  const int xw = outimg.xw;
  const int yw = outimg.yw;
  const int ndet = detpoly.size();
  const int nmask = maskpoly.size();

  // get bounds for each polygon
  std::vector<Rect> detbounds, maskbounds;
  detbounds.reserve(ndet); maskbounds.reserve(nmask);
  for(auto& p : detpoly)
    detbounds.emplace_back(p.bounds());
  for(auto& p : maskpoly)
    maskbounds.emplace_back(p.bounds());

  // list of gradients for each polygon
  std::vector< std::vector<float> > detgrads, maskgrads;
  for(auto& p : detpoly)
    detgrads.emplace_back(getGrads(p));
  for(auto& p : maskpoly)
    maskgrads.emplace_back(getGrads(p));

  // iterate over scanlines in image
  std::vector<intersect> isects;
  isects.reserve(32);

  float lastdetx=0, thisdetx=0;

  for(int y=0; y<yw; ++y)
    {
      // get at which x values the polygons cross the scan line
      const float yf = y;
      isects.clear();
      for(int i=0; i<ndet; ++i)
        if(yf >= detbounds[i].tl.y && yf <= detbounds[i].br.y)
          addCrossings(isects, yf, detpoly[i], detgrads[i], true);
      for(int i=0; i<nmask; ++i)
        if(yf >= maskbounds[i].tl.y && yf <= maskbounds[i].br.y)
          addCrossings(isects, yf, maskpoly[i], maskgrads[i], false);

      // no polygons cross this scanline, so ignore
      if(isects.empty())
        continue;

      // sort in x order
      std::sort(isects.begin(), isects.end(),
                [](const intersect& a, const intersect& b) { return a.x<b.x; });

      // iterate over scanline parts, drawing between visible lines
      bool indet=false, inmask=false;
      for(auto& isect : isects)
        {
          if(isect.det)
            {
              if(indet)
                {
                  // fill row
                  indet = false;
                  thisdetx = std::min(xw-1.f, isect.x);
                  if(!inmask)
                    doFill(y, outimg, lastdetx, thisdetx, fillval);
                }
              else
                {
                  // starting a new detector segment
                  lastdetx = std::max(0.f, isect.x);
                  indet = true;
                }
            }
          else
            {
              if(inmask)
                {
                  // ending mask - store position so we can draw to next posn
                  inmask = false;
                  lastdetx = std::max(0.f, isect.x);
                }
              else
                {
                  // starting mask - draw if necessary
                  thisdetx = std::min(xw-1.f, isect.x);
                  if(indet)
                    doFill(y, outimg, lastdetx, thisdetx, fillval);
                  inmask = true;
                }
            }
        }
    }
}
/*

int main()
{
  int xw = 512, yw= 512;
  Image<float> img(xw, yw);

  PolyVec dpolys;

  Poly p1;
  p1.add(Point(  5, 50));
  p1.add(Point( 50,150));
  p1.add(Point(150,150));
  p1.add(Point(150, 50));
  dpolys.push_back(p1);

  Poly p2;
  p2.add(Point(200, 50));
  p2.add(Point(200,150));
  p2.add(Point(350,150));
  p2.add(Point(350, 50));
  dpolys.push_back(p2);


  PolyVec mpolys;
  Poly m1;
  m1.add(Point( 20, 60));
  m1.add(Point( 20,140));
  m1.add(Point(301,140));
  m1.add(Point(401, 60));
  mpolys.push_back(m1);

  fillPoly2(dpolys, mpolys, img, 1);

  write_fits_image("out.fits", img, 0, 0, 1);
}
*/

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
