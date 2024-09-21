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

Mask::Mask(const std::string& filename)
{
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

  PolyVec polys = mask_to_polygons(maskimg, true);
  std::printf("  - found %ld polygons\n", polys.size());

  // FILE *fout = std::fopen("test.reg", "w");
  // fprintf(fout, "# Region file format: DS9 version 4.1\n");
  // fprintf(fout, "fk5\n");

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
          // fprintf(fout, "polygon(");
          // bool first = true;
          // for(auto c : skyvec)
          //   {
          //     if(!first) fprintf(fout, ",");
          //     first = false;
          //     fprintf(fout, "%.7f,%.7f", c.lon, c.lat);
          //   }
          // fprintf(fout,")\n");

          ct += skyvec.size();
          maskcoords.push_back(skyvec);
        }
    }

  // std::fclose(fout);

  std::printf("  - converted to %ld sky coordinates\n", ct);
}

PolyVec Mask::as_ccd_poly(const CoordConv& cc)
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
          poly.add(Point(ccdx, ccdy));
        }
    }

  return polys;
}
