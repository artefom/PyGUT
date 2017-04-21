import math
from . import utils
from . import errors
from . import constants
from osgeo import osr
from osgeo import gdal


class RasterGrid(object):
    def __init__(self, projection=None, geotransform=None, row_number=None, column_number=None, x_res=None, y_res=None,
                 min_x=None, min_y=None, max_x=None, max_y=None):
        """
        Class, containing definitions for raster grid
        :param projection:
        :param geotransform:
        :param row_number:
        :param column_number:
        :param x_res:
        :param y_res:
        :param min_x:
        :param min_y:
        :param max_x:
        :param max_y:
        """
        self.projection = projection
        self.x_res = x_res
        self.y_res = y_res
        self.min_x = min_x
        self.max_x = max_x
        self.min_y = min_y
        self.max_y = max_y

        # split geotransform into resolution and boundaries
        if geotransform is not None and row_number is not None and column_number is not None:
            self.min_x = geotransform[0]
            self.max_y = geotransform[3]
            self.x_res = abs(geotransform[1])
            self.y_res = abs(geotransform[5])
            self.max_x = self.min_x + column_number * self.x_res
            self.min_y = self.max_y - row_number * self.y_res

    def __str__(self):
        return "x_res:{},y_res:{},min_x:{},max_x:{},min_y:{},max_y:{}".format(self.min_x, self.max_x,
                                                                           self.min_y, self.max_y, self.x_res,
                                                                           self.y_res)

    def __repr__(self):
        return self.__str__()

    def aligned_with(self, other):
        """
        Check if one raster grid is aligned with another
        :param other:
        :return:
        """
        aligned = True
        if not self.is_comparable(other):
            aligned = False

        # use pixel sizes to calculate tollerance
        pixel_num = self.get_pixel_count(self.max_x, self.min_x, self.x_res)
        pixel_num = max(pixel_num, self.get_pixel_count(other.max_x, other.min_x, other.x_res))
        pixel_num = max(pixel_num, self.get_pixel_count(self.max_y, self.min_y, self.y_res))
        pixel_num = max(pixel_num, self.get_pixel_count(other.max_y, other.min_y, other.y_res))
        res = min(self.x_res, self.y_res)
        tolerance = 0.001 * res / pixel_num

        min_x_rounded = self.round_to_grid(self.min_x, other.min_x, self.x_res)
        # check that pixel maps exactly into same place
        if abs(min_x_rounded - self.min_x) > tolerance:
            aligned = False
        max_y_rounded = self.round_to_grid(self.max_y, other.max_y, self.y_res)
        # chack that again
        if abs(max_y_rounded - self.max_y) > tolerance:
            aligned = False

        return aligned

    def intersection(self, other):
        """
        Calculate intersection of current raster grid with another
        :param other: raster grid
        :return: new instance, intersection of raster grids
        """
        if not self.is_comparable(other):
            return None

        min_x = max(self.min_x, other.min_x)
        max_x = min(self.max_x, other.max_x)
        min_y = max(self.min_y, other.min_y)
        max_y = min(self.max_y, other.max_y)

        if min_x >= max_x or min_y >= max_y:
            msg = "Images don't intersect"
            raise errors.InvalidIntersection(msg)

        new_raster_grid = RasterGrid(projection=self.projection, x_res=self.x_res, y_res=self.y_res,
                                     min_x=min_x, min_y=min_y, max_x=max_x, max_y=max_y)

        return new_raster_grid

    def union(self, other):
        """
        Compute intersection of current raster grid with another.
        :param other: raster grid
        :return: new instance of raster grid
        """
        if not self.is_comparable(other):
            return None

        min_x = min(self.min_x, other.min_x)
        min_y = min(self.min_y, other.min_y)
        max_x = max(self.max_x, other.max_x)
        max_y = max(self.max_y, other.max_y)

        new_raster_grid = RasterGrid(projection=self.projection, x_res=self.x_res, y_res=self.y_res,
                                     min_x=min_x, min_y=min_y, max_x=max_x, max_y=max_y)
        return new_raster_grid

    def is_comparable(self, other):
        """
        Checks that boundaries of one raster grid is comparable to another.
        That is, if projection matches, and grids have same resolution
        :param other:
        :return:
        """
        comparable = True

        if not self.equal_resolution(other):
            comparable = False

        if not self.equal_projection(other):
            comparable = False
        return comparable

    def equal_resolution(self, other):
        """
        Check if one grid has same resolution as other
        :param other:
        :return:
        """
        return self.x_res == other.x_res and self.y_res == other.y_res

    def equal_projection(self, other):

        self_proj = str(self.projection) if self.projection is not None else ''
        other_proj = str(other.projection) if other.projection is not None else ''
        sr_self = osr.SpatialReference(wkt=self_proj)
        sr_other = osr.SpatialReference(wkt=other_proj)
        return bool(sr_self.IsSame(sr_other))

    def get_geotransform(self):
        """
        Get geotransform of current raster grid in gdal definitions
        :return: tuple, geotransform
        """
        geotransform = (self.min_x, self.x_res, 0.0, self.max_y, 0.0, -self.y_res)
        return geotransform

    def reproject(self, target_grid):
        """
        Reproject this raster grid into anthoer raster grid
        :param target_grid:
        :return: new instance of raster grid
        """

        # create coordinate transformation
        self_sref = osr.SpatialReference(str(self.projection))
        target_sref = osr.SpatialReference(str(target_grid.projection))
        c_trans = osr.CoordinateTransformation(self_sref, target_sref)

        # get new boundaries
        (top_left_x, top_left_y, z) = c_trans.TransformPoint(self.min_x, self.max_y)
        (bottom_left_x, bottom_left_y, z) = c_trans.TransformPoint(self.min_x, self.min_y)
        (top_right_x, top_right_y, z) = c_trans.TransformPoint(self.max_x, self.max_y)
        (bottom_right_x, bottom_right_y, z) = c_trans.TransformPoint(self.max_x, self.min_y)

        # sort them, so minimums and maximums are set correctly
        min_x = min(top_left_x, bottom_left_x)
        min_y = min(bottom_left_y, bottom_right_y)
        max_x = max(top_right_x, bottom_right_x)
        max_y = max(top_left_y, top_right_y)

        # round it up, so it lies strictly within another grid
        min_x = self.round_to_grid(min_x, target_grid.min_x, target_grid.x_res)
        min_y = self.round_to_grid(min_y, target_grid.min_y, target_grid.y_res)
        max_x = self.round_to_grid(max_x, target_grid.min_x, target_grid.x_res)
        max_y = self.round_to_grid(max_y, target_grid.min_y, target_grid.y_res)

        new_raster_grid = RasterGrid(projection=target_grid.projection, x_res=target_grid.x_res,
                                     y_res=target_grid.y_res, min_x=min_x, min_y=min_y, max_x=max_x,
                                     max_y=max_y)

        return new_raster_grid

    def get_size(self):
        """
        Get number of rows and columns of current raster grid
        :return:
        """
        rows_count = self.get_pixel_count(self.max_y, self.min_y, self.y_res)
        columns_count = self.get_pixel_count(self.max_x, self.min_x, self.x_res)
        return rows_count, columns_count

    @staticmethod
    def get_pixel_count(grid_max, grid_min, grid_res):
        """
        Get number of pixels by min, max and resolution
        :return: int, number of pixels
        """
        npix = int(round((grid_max - grid_min) / grid_res))
        return npix

    @staticmethod
    def round_to_grid(val, val_on_grid, res):
        """
        internal use only.
        rounds up val, so that is round number of pixels away from val_on_grid
        :param val:
        :param val_on_grid:
        :param res:
        :return:
        """
        diff = val - val_on_grid
        pixel_num = diff / res
        pixel_round_num = round(pixel_num)
        snapped_val = val_on_grid + pixel_round_num * res
        return snapped_val


def find_common_region(grid_list, ref_grid, footprint=constants.FP_INTERSECTION):
    """
    Finds common raster grid of rasters based on footprint enum
    :param grid_list: list of raster grids to find common region
    :param ref_grid: reference gird (return this value if footprint == FP_REFERENCE
    :param footprint: FP_REFERENCE or FP_INTERSECTION or FP_UNION
    :return:
    """
    new_grid = ref_grid
    if footprint != constants.FP_REFERENCE:
        for grid in grid_list:
            if not new_grid.aligned_with(grid):
                grid = grid.reproject(ref_grid)

            if footprint == constants.FP_INTERSECTION:
                new_grid = new_grid.intersection(grid)
            elif footprint == constants.FP_UNION:
                new_grid = new_grid.union(grid)

    return new_grid

def raster_grid_from_file(filename):
    """
    Calculate raster grid definitions based on file
    :param filename:
    :return: instance of RasterGrid
    """
    ds = gdal.Open(filename)
    geotransform = ds.GetGeoTransform()
    nrows = ds.RasterYSize
    ncols = ds.RasterXSize
    projection = ds.GetProjection()
    rastgrid = RasterGrid(projection=projection, geotransform=geotransform, row_number=nrows,
                          column_number=ncols)
    return rastgrid
