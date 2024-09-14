#ifndef GEOM_HH
#define GEOM_HH

#include <vector>

struct Point
{
  Point() : x(0), y(0) {}
  Point(float _x, float _y) : x(_x), y(_y) {}
  float x, y;
};

// tl: top-left, br: bottom-right
struct Rect
{
  Rect() {}
  Rect(Point _tl, Point _br) : tl(_tl), br(_br) {}
  Point tl, br;
};

typedef std::vector<Point> Poly;

template <typename T> T clip(T v, T minv, T maxv)
{
  return v < minv ? minv : (v > maxv ? maxv : v);
}

// clip polygons (polys must be defined the right way round)
Poly poly_clip(const Poly& spoly, const Poly& cpoly);

// area of polygon (positive: direction correct)
float poly_area(const Poly& p);

// return rectangle bounding polygon
Rect poly_bounds(const Poly& p);

#endif
