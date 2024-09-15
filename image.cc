#include <tuple>
#include "image.hh"

static constexpr int INCL=-1;
static constexpr int EXCL=-2;

// try to recursively build increasing sized rectangle starting from x0,y0
// xr, yr is the size of the already-tested rectangle
// expx/y: can we increase in size in these directions?
static std::tuple<int, int>
inner_rect(const Image<int>& img, int x0, int y0,
           int xr, int yr, bool expx, bool expy)
{
  bool okx=false, oky=false, okc=false;

  if(expx && x0+xr < int(img.xw))
    {
      okx = true;
      for(int y=y0; y<y0+yr; ++y)
        if(img(x0+xr,y) != INCL)
          {
            okx = false; break;
          }
    }
  if(expy && y0+yr < int(img.yw))
    {
      oky = true;
      for(int x=x0; x<x0+xr; ++x)
        if(img(x,y0+yr) != INCL)
          {
            oky = false; break;
          }
    }
  if(x0+xr < int(img.xw) && y0+yr < int(img.yw))
    {
      okc = img(x0+xr,y0+yr) == INCL;
    }

  if(!okx && !oky)
    return std::tuple<int,int>(xr,yr);

  if(okx && oky && okc)
    {
      xr++; yr++;
    }
  else if (okx && oky)
    {
      if( (xr+1)*yr > xr*(yr+1) )
        xr++;
      else
        yr++;
    }
  else
    {
      if(okx)
        {
          xr++; expy = false;
        }
      else if(oky)
        {
          yr++; expx = false;
        }
    }
  return inner_rect(img, x0, y0, xr, yr, expx, expy);
}


std::vector<Poly> mask_poly(const Image<int>& _img)
{
  Image<int> img(_img.xw, _img.yw);
  std::vector<Poly> out;

  // we rewrite the regions so we have space
  for(unsigned y=0; y!=img.yw; ++y)
    for(unsigned x=0; x!=img.yw; ++x)
      img(x,y) = _img(x,y) > 0 ? INCL : EXCL;
 
  //

  return out;
}
