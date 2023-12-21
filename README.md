COTAT+ algorithm for flow direction upscaling
=======
COTAT+ algorithm (Cell Outlet Tracing with an Area Threshold plus)

CONTACT
-------------

adrianorpaz@ct.ufpb.br
 
+55 83 3216 7355 

HARDWARE REQUIREMENTS
-------------
Devices that run Windows 10 for desktop editions (a 1GHz or faster processor)

SOFTWARE ESPECIFICATIONS
-------------
**PROGRAM LANGUAGE:**	Fortran 90 

**SOFTWARE REQUIREMENTS:**	Fortran 90 compiler, Windows 

**PROGRAM SIZE:**	36.2 Kb
 
AVAILABILITY
-------------
The source codes are available for download at the link: 

	https://github.com/lacmaufpb/Upscaling-1.0 

REFERENCE
============

Paz, A.R.; Collischonn, W.; Silveira, A.L.L. Improvements in large scale drainage networks derived from digital elevation models. Water Resources Research, v. 42, p. W08502, 2006. https://doi.org/10.1029/2005WR004544

VERSION
============
This algorithm version is adapted from Paz et al. (2006) for reading the high-resolution data and processing each low-resolution cell individually. It implies low RAM requirements but a long run time.

USAGE
============
The following input files are required, with the specification as described:
* **DIRHIGH (.rst;.rdc):** High-resolution flow directions. Raster file, Idrisi/TerrSet format, integer/binary. The following flow direction code must be considered:  1 (northeast), 2 (east), 4 (southeast), 8 (south), 16 (southwest), 32 (west), 64 (northwest), and 128 (north).
* **AREAHIGH (.rst;.rdc):** High-resolution flow accumulated areas. Raster file, Idrisi/TerrSet format, real/binary.
* **MASK (optional; .rst;.rdc):** Mask to indicate areas to be/not to be processed (pixels with value 1 or 0, respectively). Raster file, low resolution, Idrisi/TerrSet format, integer/binary. 
* **INPUT_UPSCALING (.txt):**  Upscaling configuration information. Ascii file.
The names of the input files are defined along the code and could be changed with caution.
The configuration of the upscaling is set in the 'input_upscaling' input file. 
The following information is defined in this file: a) high-resolution cell size (same units of the input raster); b) low-resolution cell size (same units of the input raster); c) parameter Area Threshold (AT; km2); d) parameter Minimum Upstream Flow Path (MUFP; km); e) Selection for using (1) or not (0) a mask.
For example, for flow direction upscaling from 1 m to 100 m, option a) must assume the value '1', and option b) must assume the value '100'. The definition of AT and MUFP are according to Paz et al. (2006), which recommend AT equal to the low-resolution cell area and MUFP equal to 1/5 of the low-resolution cell size. In the example, AT should be set as 0.01 (i.e., 10000 m2, but in km2), and MUFP should be set as 0.02 (i.e., 20 m, but in km).

The upscaling algorithm requires an integer number of pixels (high resolution) inside each low-resolution cell. Therefore, the ratio between the high-resolution pixel and the low-resolution cell must be an integer. Additionally, the geographic extension of all input raster must be the same. These input rasters must also follow the specifications regarding data format. The QGIS software could convert raster files from other formats to the Idrisi format required by the algorithm. The Idrisi format was adopted due to having a specific file for the data matrix (the .rst file) and another specific format for the documentation (the .rdc file), and both being very easy to read by a Fortran code. The output file follows this format and could be converted to other raster formats using the QGIS software.
The flow direction codes of input and output files could also be reclassified from or to other coding systems using the QGIS software.

The algorithm produces the following output file:
* **DIRLOW (.rst;.rdc):** Low-resolution flow directions. Raster file, Idrisi/TerrSet format, integer/binary. The same flow direction code of the DIRHIGH file is followed: 1 (northeast), 2 (east), 4 (southeast), 8 (south), 16 (southwest), 32 (west), 64 (northwest), and 128 (north).
The name of the output file is defined along the code and could be changed with caution.
