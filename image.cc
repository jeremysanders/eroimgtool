#include <limits>
#include <unistd.h>

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
    unlink(filename.c_str());

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
