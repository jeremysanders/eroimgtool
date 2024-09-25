#ifndef POLY_FILL_HH
#define POLY_FILL_HH

#include <cstdint>

#include "geom.hh"
#include "image.hh"

void fillPoly(const Poly& poly, Image<uint8_t>& outimg, uint8_t val);

void fillPoly2(const PolyVec& detpoly, const PolyVec& maskpoly,
               Image<float>& outimg, float fillval);

#endif
