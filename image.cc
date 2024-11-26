#include <cmath>
#include <limits>
#include <filesystem>
#include <stdexcept>

#include <fitsio.h>

#include "common.hh"
#include "image.hh"

// write a simple header for the image for a linear scale
static void write_header(fitsfile* ff, float xc, float yc, float pixscale)
{
  int status = 0;

  float t;
  t = xc+1;
  fits_write_key(ff, TFLOAT, "CRPIX1", &t, 0, &status);
  t = yc+1;
  fits_write_key(ff, TFLOAT, "CRPIX2", &t, 0, &status);
  t = pixscale;
  fits_write_key(ff, TFLOAT, "CDELT1", &t, 0, &status);
  fits_write_key(ff, TFLOAT, "CDELT2", &t, 0, &status);
  t = 0;
  fits_write_key(ff, TFLOAT, "CRVAL1", &t, 0, &status);
  fits_write_key(ff, TFLOAT, "CRVAL2", &t, 0, &status);

  fits_write_key(ff, TSTRING, "CUNIT1", const_cast<char*>("pix"), 0, &status);
  fits_write_key(ff, TSTRING, "CUNIT2", const_cast<char*>("pix"), 0, &status);

  check_fitsio_status(status);
}

static int get_int_bitpix(const Image<int>& img)
{
  int minv = img.arr.min();
  int maxv = img.arr.max();

  // use shortest size to hold integer
  int bitpix = LONG_IMG;
  if(minv >=0 && maxv < 256)
    bitpix = BYTE_IMG;
  else if(minv >= -128 && maxv < 128)
    bitpix = SBYTE_IMG;
  else if(minv >= -32768 && maxv < 32768)
    bitpix = SHORT_IMG;
  else if(minv >= 0 && maxv < 65536)
    bitpix = USHORT_IMG;

  return bitpix;
}


void write_fits_image(const std::string& filename,
                      const Image<int>& img,
                      float xc, float yc, float pixscale,
                      bool overwrite)
{
  if(overwrite)
    std::filesystem::remove(filename);

  int status = 0;

  fitsfile* ff;
  fits_create_file(&ff, filename.c_str(), &status);
  check_fitsio_status(status);

  long dims[] = {img.xw, img.yw};
  long fpixel[] = {1,1};
  fits_create_img(ff, get_int_bitpix(img), 2, dims, &status);
  fits_write_pix(ff, TINT, fpixel, img.xw*img.yw,
                 const_cast<int*>(&img.arr[0]),
                 &status);
  check_fitsio_status(status);
  write_header(ff, xc, yc, pixscale);

  fits_close_file(ff, &status);
  check_fitsio_status(status);
}

static std::tuple<std::valarray<int>, double> make_int_arr(const std::valarray<float> &inarr, int maxval)
{
  double arrmax = inarr.max();
  if(inarr.min() < 0)
    throw std::runtime_error("Array contains negative elements");

  // prevent division by zero
  if(arrmax == 0.)
    arrmax = 1.;

  double scale = maxval / arrmax;

  std::valarray<int> intarr(inarr.size());

  for(size_t i=0; i != intarr.size(); ++i)
    intarr[i] = int(std::round(inarr[i] * scale));

  return std::make_tuple(intarr, 1/scale);
}

void write_fits_image(const std::string& filename,
                      const Image<float>& img,
                      float xc, float yc, float pixscale,
                      bool overwrite,
                      int bitpix)
{
  if(overwrite)
    std::filesystem::remove(filename);

  int status = 0;

  fitsfile* ff;
  fits_create_file(&ff, filename.c_str(), &status);
  check_fitsio_status(status);

  long dims[] = {img.xw, img.yw};
  long fpixel[] = {1,1};

  switch(bitpix)
    {
    case -32:
      fits_create_img(ff, FLOAT_IMG, 2, dims, &status);
      fits_write_pix(ff, TFLOAT, fpixel, img.xw*img.yw,
                     const_cast<float*>(&img.arr[0]),
                     &status);
      break;
    case 8:
      {
        auto [intarr, scale] = make_int_arr(img.arr, 255);
        fits_create_img(ff, BYTE_IMG, 2, dims, &status);
        fits_write_pix(ff, TINT, fpixel, img.xw*img.yw,
                       const_cast<int*>(&intarr[0]),
                       &status);
        fits_write_key(ff, TDOUBLE, "BSCALE", &scale, "Data scaling", &status);
        double zero = 0;
        fits_write_key(ff, TDOUBLE, "BZERO", &zero, "Data offset", &status);
      }
      break;
    case 16:
      {
        auto [intarr, scale] = make_int_arr(img.arr, 32767);
        fits_create_img(ff, SHORT_IMG, 2, dims, &status);
        fits_write_pix(ff, TINT, fpixel, img.xw*img.yw,
                       const_cast<int*>(&intarr[0]),
                       &status);
        fits_write_key(ff, TDOUBLE, "BSCALE", &scale, "Data scaling", &status);
        double zero = 0;
        fits_write_key(ff, TDOUBLE, "BZERO", &zero, "Data offset", &status);
      }
      break;

    default:
      throw std::runtime_error("Invalid bitpix");
    }
      check_fitsio_status(status);
  write_header(ff, xc, yc, pixscale);

  fits_close_file(ff, &status);
  check_fitsio_status(status);
}

Image<float> read_fits_image(const std::string& filename)
{
  int status = 0;

  fitsfile* ff;
  fits_open_file(&ff, filename.c_str(), READONLY, &status);
  check_fitsio_status(status);

  int naxis;
  long naxes[2];

  fits_get_img_dim(ff, &naxis, &status);
  check_fitsio_status(status);

  if(naxis != 2)
    throw std::runtime_error("Invalid number of dimensions");

  fits_get_img_size(ff, 2, &naxes[0], &status);
  check_fitsio_status(status);

  Image<float> outimg(naxes[0], naxes[1]);

  long fpixel[] = {1,1};
  int anynul = 0;
  fits_read_pix(ff, TFLOAT, fpixel, naxes[0]*naxes[1], nullptr,
                &outimg.arr[0], &anynul, &status);
  check_fitsio_status(status);

  fits_close_file(ff, &status);
  check_fitsio_status(status);

  return outimg;
}
