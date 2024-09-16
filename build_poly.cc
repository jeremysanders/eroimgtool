#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <list>
#include <vector>

#include "image.hh"
#include "geom.hh"
#include "build_poly.hh"

// Routine to turn a mask image into a list of polygons

// The algorithm is to make a list of line segments, described by a
// starting coordinate and direction, around an initial pixel in the
// mask.

// We examine each segment. If there is a pixel included in the mask,
// but not processed on either side of the segment, then we replace
// the segment with three segments which encompass this pixel.  We
// then repeatedly examine each segment until there are no pixels to
// add around the segments.

// We then clean up the segment list to remove segments which go one
// way, then the next goes back.

// The segment list is then turned into a Polygon.

// The process is repeated until there are no more pixels included in
// the initial mask.
//////////////////////////////////////////////////////////////////////


// direction segments can be in
enum direction { RIGHT, DOWN, LEFT, UP, INVALID };

// segment specifies starting point (x,y) and direction
struct Segment
{
  Segment() {}
  Segment(int _x, int _y, direction _dir)
    : x(_x), y(_y), dir(_dir) {}
  int x, y;
  direction dir;

  // get x at end of segment
  int endx() const
  {
    switch(dir)
      {
      case LEFT:  return x-1;
      case RIGHT: return x+1;
      default:    return x;
      }
  }

  // get y at end of segment
  int endy() const
  {
    switch(dir)
      {
      case UP:   return y+1;
      case DOWN: return y-1;
      default:   return y;
      }
  }
};

// are the two segments in opposing directions?
inline bool opposing(const Segment a, const Segment b)
{
  return
    (a.dir==LEFT  && b.dir==RIGHT) ||
    (a.dir==RIGHT && b.dir==LEFT ) ||
    (a.dir==UP    && b.dir==DOWN ) ||
    (a.dir==DOWN  && b.dir==UP   );
}

static void cleanup_opposing(std::list<Segment>& segs)
{
  // clean up line segments, removing bits where we go forward,
  // then back
  for(auto s1=segs.begin(); s1 != segs.end(); )
    {
      auto s2 = std::next(s1);
      if(s2 == segs.end())
        break;
      if(opposing(*s1, *s2))
        {
          // erase the two elements
          s1 = segs.erase(s1);
          s1 = segs.erase(s1);
          if(s1 != segs.begin())
            {
              // move back one step to check if the previous was
              // also part of a pair
              s1 = std::prev(s1);
            }
        }
      else
        {
          // move to next
          s1 = s2;
        }
    }


  // move same direction segments to one end of the list, so they can
  // be coalesced when the polygon is created
  if(!segs.empty())
    {
      while(segs.front().dir == segs.back().dir)
        {
          segs.splice(segs.end(), segs, segs.begin());
        }
    }

  // also remove opposing front and back of segments
  while(!segs.empty() && opposing(segs.front(), segs.back()))
    {
      segs.pop_front();
      segs.pop_back();
    }
}

// convert series of segments to a polygon
static Poly segs_to_poly(const std::list<Segment>& segs)
{
  Poly retn;
  direction lastdir = INVALID;

  for(auto& s : segs)
    {
      // check that there's no jump in the produced segments
      assert(lastdir==INVALID || (s.x==retn.back().x && s.y==retn.back().y));

      // move segment starting point to destination
      Point endpt(s.endx(), s.endy());
      if(s.dir == lastdir)
        // go in same direction, so we replace the last point
        retn.back() = endpt;
      else
        // add new point
        retn.add(endpt);

      lastdir = s.dir;
    }

  return retn;
}

PolyVec mask_to_polygons(const Image<int>& inmask)
{
  const int xw = inmask.xw;
  const int yw = inmask.yw;

  // copy of mask, where unprocessed pixels are <0
  Image<int> mask(xw, yw);
  for(int y=0; y<yw; ++y)
    for(int x=0; x<xw; ++x)
      mask(x,y) = inmask(x,y)>0 ? -1 : 0;

  // output polygons
  std::vector<Poly> polys;

  // current list of line segments (yes, std::list is inefficient, but
  // we need to insert in middle)
  std::list<Segment> segs;

  unsigned lastidx = 0; // keep track of where the last found coordinate was
  for(int polyidx=1; ; polyidx++)
    {
      // find next pixel to start from
      {
        // scan through underlying array in 1d so we can keep track
        // of last pixel and not redo from start
        int idx=-1;
        for(unsigned i=lastidx; i != mask.size(); ++i)
          if(mask.arr[i]<0)
            {
              idx = i; break;
            }
        // nothing left
        if(idx<0)
          break;
        lastidx = idx+1;
        int x0 = idx % mask.xw;
        int y0 = idx / mask.xw;

        segs.emplace_back(x0,   y0,   UP   );
        segs.emplace_back(x0,   y0+1, RIGHT);
        segs.emplace_back(x0+1, y0+1, DOWN );
        segs.emplace_back(x0+1, y0,   LEFT );
        mask(x0,y0) = polyidx;
      }

      for(;;)
        {
          bool anyrepl=false;
          for(auto sit = segs.begin(); sit != segs.end(); )
            {
              const int x = sit->x;
              const int y = sit->y;
              Segment newseg1, newseg2, newseg3;

              // If there's a neighbouring unprocessed pixel we replace
              // the line segment with one which goes around it.
              bool repl = false; // replaced segment?
              switch(sit->dir)
                {
                case RIGHT:
                  if(y<yw && mask(x,y)<0)
                    {
                      newseg1 = Segment(x,   y,   UP);
                      newseg2 = Segment(x,   y+1, RIGHT);
                      newseg3 = Segment(x+1, y+1, DOWN);
                      mask(x,y) = polyidx;
                      repl = true;
                    }
                  else if(y>0 && mask(x,y-1)<0)
                    {
                      newseg1 = Segment(x,   y,   DOWN);
                      newseg2 = Segment(x,   y-1, RIGHT);
                      newseg3 = Segment(x+1, y-1, UP);
                      mask(x,y-1) = polyidx;
                      repl = true;
                    }
                  break;

                case DOWN:
                  if(x<xw && mask(x,y-1)<0)
                    {
                      newseg1 = Segment(x,   y,   RIGHT);
                      newseg2 = Segment(x+1, y,   DOWN);
                      newseg3 = Segment(x+1, y-1, LEFT);
                      mask(x,y-1) = polyidx;
                      repl = true;
                    }
                  else if(x>0 && mask(x-1,y-1)<0)
                    {
                      newseg1 = Segment(x,   y,   LEFT);
                      newseg2 = Segment(x-1, y,   DOWN);
                      newseg3 = Segment(x-1, y-1, RIGHT);
                      mask(x-1,y-1) = polyidx;
                      repl = true;
                    }
                  break;

                case LEFT:
                  if(y<yw && mask(x-1,y)<0)
                    {
                      newseg1 = Segment(x,   y,   UP);
                      newseg2 = Segment(x,   y+1, LEFT);
                      newseg3 = Segment(x-1, y+1, DOWN);
                      mask(x-1,y) = polyidx;
                      repl = true;
                    }
                  else if(y>0 && mask(x-1,y-1)<0)
                    {
                      newseg1 = Segment(x,   y,   DOWN);
                      newseg2 = Segment(x,   y-1, LEFT);
                      newseg3 = Segment(x-1, y-1, UP);
                      mask(x-1,y-1) = polyidx;
                      repl = true;
                    }
                  break;

                case UP:
                  if(x>0 && mask(x-1,y)<0)
                    {
                      newseg1 = Segment(x,   y,   LEFT);
                      newseg2 = Segment(x-1, y,   UP);
                      newseg3 = Segment(x-1, y+1, RIGHT);
                      mask(x-1,y) = polyidx;
                      repl = true;
                    }
                  else if(x<xw && mask(x,y)<0)
                    {
                      newseg1 = Segment(x,   y,   RIGHT);
                      newseg2 = Segment(x+1, y,   UP);
                      newseg3 = Segment(x+1, y+1, LEFT);
                      mask(x,y) = polyidx;
                      repl = true;
                    }
                  break;

                default:
                  break;
                }

              if(repl)
                {
                  // insert new segments backwards so order is correct
                  *sit = newseg3;
                  sit = segs.insert(sit, newseg2);
                  sit = segs.insert(sit, newseg1);

                  // skip forward so we're not looking always in the
                  // same direction
                  ++sit;

                  // we made a change
                  anyrepl = true;
                }

              ++sit;
            }

          // nothing left
          if(!anyrepl) break;
        }

      // remove useless opposing segments
      cleanup_opposing(segs);

      // create new polygon and add to list
      polys.emplace_back( segs_to_poly(segs) );

      segs.clear();
    }

  // for(int y=0; y<yw; ++y)
  //   {
  //     for(int x=0; x<xw; ++x)
  //       {
  //         std::printf("%d", mask(x,y));
  //       }
  //     std::printf("\n");
  //   }

  return polys;
}

/*
void test_rand()
{
  constexpr int w=10;

  Image<int> m(w,w,0);

  for(unsigned i=0; i<100; i++)
    {
      m(std::rand()%w,std::rand()%w) = 1;
    }

  for(unsigned y=0; y<m.yw; ++y)
    {
      for(unsigned x=0; x<m.xw; ++x)
        {
          std::printf("%d", m(x,y));
        }
      std::printf("\n");
    }
  std::printf("\n");

  int intarea=0;
  for(int y=0; y<w; ++y)
    for(int x=0; x<w; ++x)
      intarea += m(x,y);

  std::vector<Poly> polys(mask_to_polygons(m));
  float totarea = 0;
  for(auto poly : polys)
    {
      totarea += poly_area(poly);
    }
  std::printf("num=%ld, total area=%g, intarea=%d\n", polys.size(), totarea, intarea);
  if(totarea-intarea != 0.f)
    std::printf("BAD!!!!\n");
}

int main()
{
  for(int i=0; i<100; ++i)
    test_rand();
}
*/

/*
int main()
{
  Image<int> m(6,5);
  m(0,0)=0; m(1,0)=1; m(2,0)=1; m(3,0)=0; m(4,0)=0; m(5,0)=1;
  m(0,1)=1; m(1,1)=1; m(2,1)=1; m(3,1)=1; m(4,1)=0; m(5,1)=1;
  m(0,2)=0; m(1,2)=1; m(2,2)=0; m(3,2)=1; m(4,2)=0; m(5,2)=0;
  m(0,3)=0; m(1,3)=1; m(2,3)=1; m(3,3)=1; m(4,3)=0; m(5,3)=0;
  m(0,4)=0; m(1,4)=0; m(2,4)=0; m(3,4)=0; m(4,4)=1; m(5,4)=0;

  std::vector<Poly> polys(mask_to_polygons(m));

  for(auto poly : polys)
    {
      std::printf("size=%g\n", poly_area(poly));
      for(auto pt : poly)
        {
          std::printf("%g,%g\n", pt.y, pt.x);
        }
      std::printf("\n");
    }
}
*/
