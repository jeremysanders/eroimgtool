#include <wcslib/wcslib.h>
// yuck
#undef PI

#include <cstdlib>
#include <cstring>
#include <stdexcept>

#include <fitsio.h>

#include "build_poly.hh"
#include "common.hh"
#include "mask.hh"

namespace
{

  // little WCS wrapper
  class WCS
  {
  public:
    WCS(fitsfile *ff);
    ~WCS();
    CoordVec pix2sky(const CoordVec& pvec);

  private:
    int nwcs;
    struct wcsprm *wcs;
  };

  WCS::WCS(fitsfile *ff)
  {
    int status = 0;
    int nkeyrec;
    char* hdrstr;
    fits_hdr2str(ff, 0, NULL, 0, &hdrstr, &nkeyrec, &status);
    check_fitsio_status(status);

    int nreject;
    status = wcspih(hdrstr, nkeyrec, WCSHDR_all, 2, &nreject, &nwcs, &wcs);
    std::free(hdrstr);

    if(status != 0)
      {
        std::string err = std::string("WCS wcspih ERROR ") +
          std::to_string(status) + ": " + wcshdr_errmsg[status];
        throw std::runtime_error(err);
      }

    status = wcsset(wcs);
    if (status != 0)
      {
        std::string err = std::string("WCS wcsset ERROR ") +
          std::to_string(status) + ": " + wcshdr_errmsg[status];
        throw std::runtime_error(err);
      }
  }

  CoordVec WCS::pix2sky(const CoordVec& pixcrd)
  {
    int ncoord = pixcrd.size();
    if(ncoord<=0)
      return CoordVec();

    CoordVec world, imgcrd;
    world.resize(ncoord);
    imgcrd.resize(ncoord);

    std::vector<double> phi, theta;
    phi.resize(ncoord);
    theta.resize(ncoord);
    std::vector<int> stat;
    stat.resize(ncoord);

    int status = wcsp2s(wcs, ncoord, 2, &pixcrd[0].lon, &imgcrd[0].lon,
                        &phi[0], &theta[0], &world[0].lon, &stat[0]);
    if( status != 0 )
      {
        std::string err = std::string("WCS wcsp2s ERROR ") +
          std::to_string(status) + ": " + wcshdr_errmsg[status];
        throw std::runtime_error(err);
      }

    return world;
  }

  WCS::~WCS()
  {
    wcsvfree(&nwcs, &wcs);
  }

} // namespace

Mask::Mask()
{
}

Mask::Mask(const std::string& filename, bool simplify)
{
  if(filename.empty())
    {
      // ignore any empty mask
      return;
    }

  int status = 0;
  fitsfile *ff;

  std::printf("Opening mask %s\n", filename.c_str());
  fits_open_file(&ff, filename.c_str(), READONLY, &status);
  check_fitsio_status(status);

  int naxis;
  fits_get_img_dim(ff, &naxis, &status);
  check_fitsio_status(status);
  if(naxis != 2)
    {
      throw std::runtime_error("invalid number of dimensions in mask " + filename);
    }

  long axes[2];
  fits_get_img_size(ff, 2, axes, &status);
  check_fitsio_status(status);
  std::printf("  - dimensions %ld x %ld\n", axes[0], axes[1]);

  Image<int> maskimg(axes[0], axes[1]);
  long fpixel[2] = {1,1};
  fits_read_pix(ff, TINT, fpixel, axes[0]*axes[1], 0, &maskimg.arr[0],
                0, &status);
  check_fitsio_status(status);

  WCS wcs(ff);

  fits_close_file(ff, &status);
  check_fitsio_status(status);

  PolyVec polys = mask_to_polygons(maskimg, true, !simplify);
  std::printf("  - found %ld polygons\n", polys.size());

  size_t ct = 0;
  for(auto &poly : polys)
    {
      CoordVec pixvec;
      pixvec.reserve(poly.size());
      for(auto pt : poly.pts)
        pixvec.emplace_back(pt.x+0.5, pt.y+0.5);

      CoordVec skyvec = wcs.pix2sky(pixvec);
      if(skyvec.size() != 0)
        {
          ct += skyvec.size();
          maskcoords.push_back(skyvec);
        }
    }

  std::printf("  - converted to %ld sky coordinates\n", ct);


  if(simplify)
    {
      simplifyPolys();

      ct = 0;
      for(auto &coords : maskcoords)
        ct += coords.size();

      std::printf("  - simplified to %ld coordinates\n", ct);
    }

}

void Mask::setMaskPts(const std::vector<double>& pts)
{
  if(pts.size() % 3 != 0)
    throw std::runtime_error("List of parameters for masked points must "
                             "be multiple of three (ra,dec,rad_pix)");
  mask_pts = pts;
  for(unsigned i=0; i<pts.size(); i+=3)
    {
      std::printf("  - masking source (%g,%g) to radius %g pix\n",
                  pts[i], pts[i+1], pts[i+2]);
    }
}

void Mask::simplifyPolys()
{
  for(auto& cv : maskcoords)
    {
      if(cv.size() < 6)
        continue;

      CoordVec out;
      for(size_t i=0; i < cv.size(); i+=2)
        {
          Coord c1 = cv[i];
          if(i == cv.size()-1)
            out.push_back(c1);
          else
            {
              Coord c2 = cv[i+1];
              out.emplace_back(0.5*(c1.lon+c2.lon), 0.5*(c1.lat+c2.lat));
            }
        }
      cv = out;
    }
}

void Mask::writeRegion(const std::string& filename) const
{
  // no error checking! debugging only
  FILE *fout = std::fopen(filename.c_str(), "w");
  fprintf(fout, "# Region file format: DS9 version 4.1\n");
  fprintf(fout, "fk5\n");

  for(auto& cv : maskcoords)
    {
      fprintf(fout, "polygon(");
      bool first = true;
      for(auto c : cv)
        {
          if(!first) fprintf(fout, ",");
          first = false;
          fprintf(fout, "%.7f,%.7f", c.lon, c.lat);
        }
      fprintf(fout,")\n");
    }

  std::fclose(fout);
}

PolyVec Mask::as_ccd_poly(const CoordConv& cc) const
{
  PolyVec polys;

  for(auto& cv : maskcoords)
    {
      polys.emplace_back();
      Poly& poly = polys.back();
      poly.pts.reserve(cv.size());
      for(auto& coord : cv)
        {
          auto [ccdx, ccdy] = cc.radec2ccd(coord.lon, coord.lat);
          poly.pts.emplace_back(ccdx, ccdy);
        }
    }

  // add masks for sources, if requested
  for(unsigned i=0; i<mask_pts.size(); i+=3)
    {
      double ra = mask_pts[i];
      double dec = mask_pts[i+1];
      double rad = mask_pts[i+2];

      const int npts = 32; // number of points in "circular" polygon
      auto [ccdx, ccdy] = cc.radec2ccd(ra, dec);
      polys.emplace_back();
      Poly& poly = polys.back();
      poly.pts.reserve(npts);
      for(int i=0; i<npts; ++i)
        {
          double theta = (2*PI/npts) * (i+0.11);
          poly.add(Point(ccdx + rad*std::cos(theta),
                         ccdy + rad*std::sin(theta)));
        }
    }

  return polys;
}
