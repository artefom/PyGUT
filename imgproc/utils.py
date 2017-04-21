import numpy
from osgeo import gdalconst

from . import errors

numpy_gdal_type_mapping = [
    (numpy.uint8, gdalconst.GDT_Byte),
    (numpy.bool, gdalconst.GDT_Byte),
    (numpy.int16, gdalconst.GDT_Int16),
    (numpy.uint16, gdalconst.GDT_UInt16),
    (numpy.int32, gdalconst.GDT_Int32),
    (numpy.uint32, gdalconst.GDT_UInt32),
    (numpy.single, gdalconst.GDT_Float32),
    (numpy.float, gdalconst.GDT_Float64),
    (numpy.complex64, 8),
]


def gdal_type_to_numpy_type(gdaltype):
    for (numpy_type, test_gdal_type) in numpy_gdal_type_mapping:
        if test_gdal_type == gdaltype:
            return numpy_type
    raise errors.InvalidType("Unknown GDAL datatype: {}".format(gdaltype))


def numpy_type_to_gdal_type(numpytype):
    for (test_numpy_type, gdaltype) in numpy_gdal_type_mapping:
        if test_numpy_type == numpytype:
            return gdaltype
    raise errors.InvalidType("Unknown numpy datatype: {}".format(numpytype))

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return "Point({},{})".format(self.x, self.y)


def world_to_pixel(transform, geox, geoy):
    x = (transform[0] * transform[5] -
         transform[2] * transform[3] + transform[2] * geoy -
         transform[5] * geox) / (transform[2] * transform[4] - transform[1] * transform[5])

    y = (transform[1] * transform[3] - transform[0] * transform[4] -
         transform[1] * geoy + transform[4] * geox) / (transform[2] * transform[4] - transform[1] * transform[5])

    return Point(x, y)


def pixel_to_world(transform, x, y):
    geox = transform[0] + transform[1] * x + transform[2] * y
    geoy = transform[3] + transform[4] * x + transform[5] * y

    return Point(geox, geoy)