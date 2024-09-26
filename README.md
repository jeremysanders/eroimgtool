# eROImageTool

This is a tool to construct image and non-vignetted exposure maps from
eROSITA event files for calibration purposes, such as PSF
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


    Make eROSITA unvignetted detector exposure maps
    Usage: build/eroimgtool [OPTIONS] mode event image

    Positionals:
      mode ENUM:value in {expos->1,image->0} OR {1,0} REQUIRED
                                  Program mode
      event TEXT:FILE REQUIRED    Event filename
      image TEXT REQUIRED         Output image filename

    Options:
      -h,--help                   Print this help message and exit
      --proj ENUM:value in {det->2,fov->0,fov_sky->1} OR {2,0,1} [0]
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
      --threads UINT [1]          Number of threads
