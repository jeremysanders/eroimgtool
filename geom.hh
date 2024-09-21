#ifndef GEOM_HH
#define GEOM_HH

#include <cstddef>
#include <vector>

struct Point
{
  Point() : x(0), y(0) {}
  Point(float _x, float _y) : x(_x), y(_y) {}

  void operator+=(Point o) { x += o.x; y += o.y; }
  void operator-=(Point o) { x -= o.x; y -= o.y; }
  void operator*=(float v) { x *= v; y *= v; }
  void operator/=(float v) { x /= v; y /= v; }
  Point operator+(Point o) const { return Point(x+o.x,y+o.y); }
  Point operator-(Point o) const { return Point(x-o.x,y-o.y); }
  Point operator*(float v) const { return Point(x*v, y*v); }
  Point operator/(float v) const { return Point(x/v, y/v); }

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
  bool is_inside(Point pt) const;

  void operator+=(Point pt) { for(auto& p : pts) p += pt; }
  void operator-=(Point pt) { for(auto& p : pts) p -= pt; }
  void operator*=(float v) { for(auto& p : pts) p *= v; }
  void operator/=(float v) { for(auto& p : pts) p /= v; }
  Poly operator+(Point pt) const
  {
    Poly q; q.pts.reserve(size()); for(auto p : pts) q.add(p+pt); return q;
  }
  Poly operator-(Point pt) const
  {
    Poly q; q.pts.reserve(size()); for(auto p : pts) q.add(p-pt); return q;
  }
  Poly operator*(float v) const
  {
    Poly q; q.pts.reserve(size()); for(auto p : pts) q.add(p*v); return q;
  }
  Poly operator/(float v) const
  {
    Poly q; q.pts.reserve(size()); for(auto p : pts) q.add(p/v); return q;
  }

  // bounding box of polygon
  Rect bounds() const;

  // area (-ve if anticlockwise)
  float area() const;

  // rotate polygon
  void rotate(float theta);

  // points stored as simple vector
  std::vector<Point> pts;
};

typedef std::vector<Poly> PolyVec;

// clip polygons (polys must be defined the right way round)
// opoly is overwritten (not returned, so we don't have to reallocate)
void poly_clip(const Poly& spoly, const Poly& cpoly, Poly& opoly);

#endif
