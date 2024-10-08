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
  Point operator-() const { return Point(-x, -y); }

  float x, y;
};

// tl: top-left, br: bottom-right
// top-left has min coordinates, while bot-right has max coordinates
struct Rect
{
  Rect() {}
  Rect(Point _tl, Point _br) : tl(_tl), br(_br) {}
  bool overlap(const Rect& o) const
  {
    return (tl.x < o.br.x) && (br.x > o.tl.x) &&
      (tl.y < o.br.y) && (br.y > o.tl.y);
  }
  bool inside(const Point& pt) const
  {
    return (pt.x >= tl.x) && (pt.x <= br.x) && (pt.y >= tl.y) & (pt.y <= br.y);
  }

  Point tl, br;
};


// for defining a 2d rotation
struct RotationMatrix
{
  RotationMatrix() : m00(1), m01(0), m10(0), m11(1) {};
  RotationMatrix(float v00, float v01, float v10, float v11)
    : m00(v00), m01(v01), m10(v10), m11(v11) {};
  Point rotate(Point pt) const { return Point(pt.x*m00+pt.y*m01,
                                              pt.x*m10+pt.y*m11); }
  Point rotaterev(Point pt) const { return Point( pt.x*m00-pt.y*m01,
                                                 -pt.x*m10+pt.y*m11); }

  float m00, m01, m10, m11;
};

struct Poly
{
  Poly() {}
  Poly(size_t n) : pts(n) {}
  Poly(const std::vector<Point> _pts) : pts(_pts) {}

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

public:
  // points stored as simple vector
  std::vector<Point> pts;
};

typedef std::vector<Poly> PolyVec;

// is point inside any of the polygons?
inline bool is_inside(const PolyVec& pv, const Point& pt)
{
  for(auto& poly : pv)
    if( poly.is_inside(pt) )
      return true;
  return false;
}

// clip polygons (polys must be defined the right way round)
// opoly is overwritten (not returned, so we don't have to reallocate)
void poly_clip(const Poly& spoly, const Poly& cpoly, Poly& opoly);

// move to origin origrot, apply rotation matrix, then scale output, then shift to origimg
void applyShiftRotationScaleShift(PolyVec& polys, const RotationMatrix& mat, Point origrot,
                                  float scale, Point origimg);

#endif
