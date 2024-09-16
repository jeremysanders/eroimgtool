#ifndef GEOM_HH
#define GEOM_HH

#include <vector>

struct Point
{
  Point() : x(0), y(0) {}
  Point(float _x, float _y) : x(_x), y(_y) {}

  void operator+=(Point o) { x += o.x; y += o.y; }
  void operator-=(Point o) { x -= o.x; y -= o.y; }
  void operator*=(float v) { x *= v; y *= v; }
  Point operator+(Point o) const { return Point(x+o.x,y+o.y); }
  Point operator-(Point o) const { return Point(x-o.x,y-o.y); }
  Point operator*(float v) const { return Point(x*v, y*v); }

  float x, y;
};

// tl: top-left, br: bottom-right
struct Rect
{
  Rect() {}
  Rect(Point _tl, Point _br) : tl(_tl), br(_br) {}
  Point tl, br;
};

struct Poly
{
  Poly() {}
  Poly(const std::vector<Point> _pts)
    : pts(_pts)
  {}
  Point& operator[](size_t idx) { return pts[idx]; }
  Point operator[](size_t idx) const { return pts[idx]; }
  size_t size() const { return pts.size(); }
  void add(Point pt) { pts.push_back(pt); }
  void clear() { pts.clear(); }
  Point& back() { return pts.back(); }
  Point& front() { return pts.front(); }
  bool empty() const { return pts.empty(); }

  // bounding box of polygon
  Rect bounds() const;

  // area (-ve if anticlockwise)
  float area() const;

  void operator+=(Point pt) { for(auto p : pts) p += pt; }
  void operator-=(Point pt) { for(auto p : pts) p -= pt; }
  void operator*=(float v) { for(auto p : pts) p *= v; }

  std::vector<Point> pts;
};

typedef std::vector<Poly> PolyVec;

template <typename T> T clip(T v, T minv, T maxv)
{
  return v < minv ? minv : (v > maxv ? maxv : v);
}

// clip polygons (polys must be defined the right way round)
Poly poly_clip(const Poly& spoly, const Poly& cpoly);

#endif
