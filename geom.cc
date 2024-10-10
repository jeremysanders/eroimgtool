#include "geom.hh"
#include <algorithm>
#include <cmath>
#include <limits>

static inline bool is_inside(Point p1, Point p2, Point q)
{
  float R = (p2.x-p1.x)*(q.y-p1.y) - (p2.y-p1.y)*(q.x-p1.x);
  return R <= 0;
}

static inline Point compute_intersection(Point p1, Point p2, Point p3, Point p4)
{
  float x, y;

  if( std::abs(p2.x-p1.x) < 1e-5f )
    {
      x = p1.x;
      float m2 = (p4.y-p3.y) / (p4.x-p3.x);
      float b2 = p3.y - m2*p3.x;
      y = m2*x + b2;
    }
  else if( std::abs(p4.x-p3.x) < 1e-5f )
    {
      x = p3.x;
      float m1 = (p2.y-p1.y) / (p2.x-p1.x);
      float b1 = p1.y - m1*p1.x;
      y = m1*x + b1;
    }
  else
    {
      float m1 = (p2.y-p1.y) / (p2.x-p1.x);
      float b1 = p1.y - m1*p1.x;
      float m2 = (p4.y-p3.y) / (p4.x-p3.x);
      float b2 = p3.y - m2*p3.x;
      x = (b2-b1) / (m1-m2);
      y = m1*x + b1;
    }

  return Point(x,y);
}

// Sutherland-Hodgman algorithm
// see https://github.com/mhdadk/sutherland-hodgman/blob/main/SH.py
// points need to be in correct order (anticlockwise?) or overlap is zero

void poly_clip(const Poly& spoly, const Poly& cpoly, Poly& opoly)
{
  opoly = spoly;
  Poly npoly;
  npoly.pts.reserve(std::max(spoly.size(), cpoly.size())*2);

  for(size_t i=0; i != cpoly.size(); ++i)
    {
      npoly = opoly;
      opoly.clear();

      Point cedge1 = cpoly[i==0 ? cpoly.size()-1 : i-1];
      Point cedge2 = cpoly[i];

      for(size_t j=0; j != npoly.size(); ++j)
        {
          Point sedge1 = npoly[j==0 ? npoly.size()-1 : j-1];
          Point sedge2 = npoly[j];
          if(is_inside(cedge1, cedge2, sedge2))
            {
              if(!is_inside(cedge1, cedge2, sedge1))
                {
                  Point inter = compute_intersection(sedge1, sedge2, cedge1, cedge2);
                  opoly.add(inter);
                }
              opoly.add(sedge2);
            }
          else if(is_inside(cedge1, cedge2, sedge1))
            {
              Point inter = compute_intersection(sedge1, sedge2, cedge1, cedge2);
              opoly.add(inter);
            }

        } // npoly

    } // cpoly
}

float Poly::area() const
{
  if(size()<3) return 0;

  float a = 0;
  size_t j = size()-1;

  for(size_t i=0; i != size(); ++i)
    {
      a += (pts[j].x+pts[i].x) * (pts[j].y-pts[i].y);
      j = i;
    }

  return a*0.5f;
}

Rect Poly::bounds() const
{
  float minx=std::numeric_limits<float>::infinity();
  float miny=std::numeric_limits<float>::infinity();
  float maxx=-std::numeric_limits<float>::infinity();
  float maxy=-std::numeric_limits<float>::infinity();

  for(auto pt : pts)
    {
      minx = std::min(minx, pt.x);
      maxx = std::max(maxx, pt.x);
      miny = std::min(miny, pt.y);
      maxy = std::max(maxy, pt.y);
    }

  return Rect(Point(minx,miny), Point(maxx,maxy));
}

void Poly::rotate(float theta)
{
  const float s = std::sin(theta);
  const float c = std::cos(theta);
  for(Point& p : pts)
    {
      float tx=p.x, ty=p.y;
      p.x = tx*c - ty*s;
      p.y = tx*s + ty*c;
    }
}

bool Poly::is_inside(Point pt) const
{
  size_t num = size();
  if(num < 3)
    return false;

  // first check bounds
  Rect b = bounds();
  if(pt.x < b.tl.x || pt.y < b.tl.y || pt.x > b.br.x || pt.y > b.br.y)
    return false;

  // use ray casting algorithm.
  // odd times horizontal ray crosses poly edge says whether inside
  unsigned count = 0;
  for(size_t i=0; i != num; ++i)
    {
      Point p1 = pts[i];
      Point p2 = pts[(i+1) % num];

      if( pt.y >= std::min(p1.y, p2.y) && pt.y <= std::max(p1.y, p2.y) )
        {
          if( pt.x < std::min(p1.x, p2.x) )
            // must cross
            ++count;
          else if( pt.x > std::max(p1.x, p2.x) )
            // can't cross
            ;
          else if( std::abs(p1.y - p2.y) > 1e-6f )
            {
              // calculate line x at y=pt.y
              float grad = (p2.x - p1.x) / (p2.y - p1.y);
              float lx = p1.x + grad * (pt.y-p1.y);
              // does +x extending line from pt cross this?
              if( lx > pt.x )
                ++count;
            }
        }
    }

  return count % 2 != 0;
}

void applyShiftRotationShift(PolyVec& polys, const Matrix2& mat,
                             Point origrot, Point origimg)
{
  for(auto& poly : polys)
    for(size_t i=0; i!=poly.size(); ++i)
      {
        float dx = poly[i].x - origrot.x;
        float dy = poly[i].y - origrot.y;
        poly[i].x = dx*mat.m00 + dy*mat.m01 + origimg.x;
        poly[i].y = dx*mat.m10 + dy*mat.m11 + origimg.y;
      }
}


/*
#include <iostream>
int main()
{
  Poly a;
  a.add(Point(1,2));
  a.add(Point(2,2));
  a.add(Point(2,1));
  a.add(Point(1,1));

  Poly b;
  b.add(Point(1,3));
  b.add(Point(3,3));
  b.add(Point(3,1));
  b.add(Point(1,1));

  Poly c = poly_clip(a, b);
  std::cout << poly_area(a) << '\n';
  std::cout << poly_area(b) << '\n';
  std::cout << poly_area(c) << '\n';
  for(auto p : c)
    {
      std::cout << p.x << ' ' << p.y << '\n';
    }
}

*/
