#ifndef BUILD_POLY_HH
#define BUILD_POLY_HH

#include "image.hh"

// invert: find 0 regions in map, not 1
// merge: join redundant line segments in same direction
PolyVec mask_to_polygons(const Image<int>& img, bool invert=false, bool merge=true);

#endif
