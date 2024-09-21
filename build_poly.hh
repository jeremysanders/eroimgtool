#ifndef BUILD_POLY_HH
#define BUILD_POLY_HH

#include "image.hh"

PolyVec mask_to_polygons(const Image<int>& img, bool invert=false);

#endif
