# Artyom Fomenko BSE142

import sys

FP_INTERSECTION = 0 # resulting footprint is intersection of inputs
FP_UNION = 1 # resulting footprint is union of inputs
FP_REFERENCE = 2 # resulting footprint is determined based on reference

# other options are:
#  'near', 'bilinear', 'cubic', 'cubicspline', 'lanczos'.
#   see gdalwarp for more info
DEFAULT_RESAMPLEMETHOD = "near"

DEFAULT_FOOTPRINT = FP_INTERSECTION
DEFAULT_WINDOW_SIZE_X = 256 # default block size
DEFAULT_WINDOW_SIZE_Y = 256 # default block size
DEFAULT_PADDING = 0
DEFAULT_LOGSTREAM = sys.stdout

# default driver for output drivers
# see gdal documentation for other drivers
DEFAULT_DRIVER_NAME = "GTiff"

default_driver_parameters = dict()
default_driver_parameters['GTiff'] = ['TILED=YES', 'COMPRESS=LZW', 'INTERLEAVE=BAND', 'BIGTIFF=IF_SAFER']

GDALWARP = 'gdalwarp'
