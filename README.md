# eROImageTool

This is a tool to construct images, pseudo event files and
non-vignetted exposure maps from eROSITA event files for calibration
purposes, such as PSF modelling. The output image is by default in
relative detector coordinates around a source position. The user can
specify a mask image in sky coordinates to remove nearby sources.

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
      mode ENUM:value in {event->2,expos->1,image->0} OR {2,1,0} REQUIRED
                                  Program mode
      event TEXT:FILE REQUIRED    Event filename
      image TEXT REQUIRED         Output image filename

    Options:
      -h,--help                   Print this help message and exit
      --sources [FLOAT,FLOAT] REQUIRED
                                  List of RA,Dec for sources
      --proj ENUM:value in {box->5,det->2,fov->0,fov_sky->1,radial->3,radial_sym->4} OR {5,2,0,1,3,4} [0]
                                  Projection mode
      --proj-args FLOAT ...       List of arguments for projection
      --tm INT:INT in [1 - 7] [1]
                                  TM number
      --pixsize FLOAT [1]         Pixel size (detector pixels)
      --mask TEXT:FILE            Input mask filename
      --mask-pts [FLOAT,FLOAT,FLOAT] ...
                                  Extra masks (list ra,dec,rad_pix)
      --detmap                    Add CALDB DETMAP mask
      --gti TEXT:FILE             Additional GTI file to merge
      --xw UINT [512]             X output image size
      --yw UINT [512]             Y output image size
      --pi-min FLOAT [300]        Minimum PI value (image/event mode)
      --pi-max FLOAT [2300]       Maximum PI value (image/event mode)
      --delta-t FLOAT [0.01]      Time step (s)
      --threads UINT [1]          Number of threads
      --bitpix INT [-32]          How many bitpix to use for output exposure maps

## Modes

  * `image`: Write an output image file containing the projected number of counts in each pixel
  * `expos`: Write an output exposure map image containing the non-vignetted exposure time in each pixel
  * `event`: Write transformed events to a FITS table. The table (HDU name EROEVT) has three columns DX, DY and PI. DX and DY are the transformed coordinates relative to the source in detector pixels. PI is taken from the input event file.

## Projection modes

  * `full`: Use all photons and time periods. The source is at centre of image, with the output in relative detector coordinates. You will also need the `--detmap` option to match standard eROSITA evtool/expmap behaviour.
  * `fov`: Source is in the standard circular detector mask. The source is located at the centre of the image, with the output in relative detector coordinates.
  * `fov_sky`: Source is in the standard circular detector mask. The source is located at the centre of the image, with the output rotated into sky coordinates.
  * `det`: Detector coordinates, with the image centred on the centre of the detector.
  * `radial`: Region within the radial range (Rin<=r<Rout) given by `--proj-args Rin Rout`.
  * `radial_sym`: Region within the radial range (Rin<=r<Rout) given by `--proj-args Rin Rout`, rotated by angle of the source on the detector, to make a symmetric PSF for a radial region.
  * `box`: Region within detector coordinates (x1<=x<x2, y1<=y<y2), given by `--proj-args x1 y1 x2 y2`.

## Example command line

    eroimgtool --tm=2 \
      --sources 255.7054905,-48.7900230 \
      --pixsize=0.25 \
      --mask=em01_258138_020_ML00001_004_c946/030_mask_final.fits.gz \
      --threads=4 \
      --proj=radial_sym --proj-args 90 180 \
      image \
      em01_258138_020_ML00001_004_c946/evt.fits.gz \
      image_out.fits
