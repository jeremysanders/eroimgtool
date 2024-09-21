Requirements
 - cfitsio
 - wcslib

Usage

eroimgtool MODE EVENT MASK PIXSIZE NPIX OUTIMG TMNR RA DEC PIMIN PIMAX PSFMODE [...]

Parameters

 * MODE - "image" or "expos"

 * EVENT - event file

 * MASK - mask image

 * PIXSIZE - pixel size in detector pixels

 * NPIX - size of output image

 * TMNR - 1..7

 * RA, DEC - coordinate of source

 * PIMIN, PIMAX - range of PI values in image mode

PSFMODE is

 * av_fov - average when src within standard FoV
 * radial minrad maxrad - average when src within radial range (pix)
 * box minx maxx miny maxy - average when src insize ccd coordinate box (pix)
