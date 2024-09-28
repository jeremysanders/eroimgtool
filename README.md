# eROImageTool

This is a tool to construct images and non-vignetted exposure maps
from eROSITA event files for calibration purposes, such as PSF
modelling. The output image is by default in relative detector
coordinates around a source position. The user can specify a mask
image in sky coordinates to remove nearby sources.

## Requirements
 - C++17 compiler (e.g. gcc)
 - [CFITSIO](https://heasarc.gsfc.nasa.gov/docs/software/fitsio/)
 - [WCSLIB](https://www.atnf.csiro.au/people/mcalabre/WCS/)

## Included external components
 - [CLI11](https://github.com/CLIUtils/CLI11)

## Building
 - Use `make` to build from `Makefile`
 - You may need to add directories to find cfitsio/wcslib
 - Output executable is `build/eroimgtool`

## Current parameters

    Make eROSITA unvignetted detector exposure maps and images
    Usage: build/eroimgtool [OPTIONS] mode event image

    Positionals:
      mode ENUM:value in {expos->1,image->0} OR {1,0} REQUIRED
                                  Program mode
      event TEXT:FILE REQUIRED    Event filename
      image TEXT REQUIRED         Output image filename

    Options:
      -h,--help                   Print this help message and exit
      --proj ENUM:value in {box->5,det->2,fov->0,fov_sky->1,radial->3,radial_sym->4} OR {5,2,0,1,3,4} [0] 
                                  Projection mode
      --proj-args FLOAT ...       List of arguments for projection
      --tm INT:INT in [1 - 7] [1]
                                  TM number
      --ra FLOAT REQUIRED         RA of source (deg)
      --dec FLOAT REQUIRED        Dec of source (deg)
      --pixsize FLOAT [1]         Pixel size (detector pixels)
      --mask TEXT:FILE            Input mask filename
      --xw UINT [512]             X output image size
      --yw UINT [512]             Y output image size
      --pi-min FLOAT [300]        Minimum PI value (image mode)
      --pi-max FLOAT [2300]       Maximum PI value (image mode)
      --delta-t FLOAT [0.01]      Time step (s)
      --mask-src-rad FLOAT [0]    Source mask radius (pix)
      --threads UINT [1]          Number of threads

## Projection modes

  * `fov`: Region within the standard circular mask, in detector coordinates. The source is located at the centre of the image.
  * `fov_sky`: Region within the standard circular mask, rotated into sky coordinates. The source is located at the centre of the image.
  * `det`: Detector coordinates, with the image centred on the centre of the detector.
  * `radial`: Region within the radial range (Rin<=r<Rout) given by `--proj-args Rin Rout`.
  * `radial_sym`: Region within the radial range (Rin<=r<Rout) given by `--proj-args Rin Rout`, rotated by angle of the source on the detector, to make a symmetric PSF for a radial region.
  * `box`: Region within detector coordinates (x1<=x<x2, y1<=y<y2), given by `--proj-args x1 y1 x2 y2`.

## Example command line

    eroimgtool --tm=2 \
      --ra=255.7054905 --dec=-48.7900230 \
      --pixsize=0.25 \
      --mask=em01_258138_020_ML00001_004_c946/030_mask_final.fits.gz \
      --threads=4 \
      --proj=radial_sym --proj-args 90 180 \
      image \
      em01_258138_020_ML00001_004_c946/evt.fits.gz \
      image_out.fits
