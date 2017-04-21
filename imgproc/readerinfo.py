import math

import numpy

from . import utils


class ReaderInfo(object):
    """
    Info object that contains all usefull information about processing process.
    Passed to user function.

    fields:

    block_size_x
    block_size_y
    padding
    logstream
    reference_grid raster grid that represents current projection, and transformation
    x_size size of reference_grid
    y_size size of reference_grid
    x_total_blocks
    y_total_blocks
    block_width
    block_height
    block_lookup lookup table that is used to derive dataset of each block
    """
    def __init__(self, reference_grid,
                 block_size_x, block_size_y,
                 padding, logstream):

        self.block_size_x = block_size_x
        self.block_size_y = block_size_y
        self.padding = padding

        self.logstream = logstream

        self.reference_grid = reference_grid

        self.x_size = int(round((self.reference_grid.max_x -
                                 self.reference_grid.min_x) / self.reference_grid.x_res))
        self.y_size = int(round((self.reference_grid.max_y -
                                 self.reference_grid.min_y) / self.reference_grid.y_res))

        self.x_total_blocks = int(math.ceil(float(self.x_size) / self.block_size_x))
        self.y_total_blocks = int(math.ceil(float(self.y_size) / self.block_size_y))

        self.block_padded_width = None
        self.block_padded_height = None

        self.block_top_left = None
        self.block_bottom_right = None

        self.x_block = None
        self.y_block = None

        self.block_lookup = {}

    def set_block_dataset(self, block, dataset, filename):
        """
        Assing gdal dataset to specific block
        :param block: numpy array
        :param dataset: gdal dataset
        :param filename: path to dataset
        :return:
        """
        self.block_lookup[id(block)] = (dataset, filename)

    def get_block_size(self):
        """
        Get size of block (without padding)
        :return: self.block_size_x, self.block_size_y
        """
        return self.block_size_x, self.block_size_y

    def get_padding_size(self):
        """
        Get size of padding
        :return: self.padding
        """
        return self.padding

    def get_total_size(self):
        """
        Get size of processed area
        :return: self.x_size, self.y_size
        """
        return self.x_size, self.y_size

    def get_geotransform(self):
        """
        Get geotransform in gdal definitions
        :return: self.reference_grid.get_geotransform()
        """
        return self.reference_grid.get_geotransform()

    def get_projection(self):
        """
        Get projection of current working region
        :return: self.reference_grid.projection
        """
        return self.reference_grid.projection

    def get_total_blocks(self):
        """
        Get total number of blocks
        :return: x_total_blocks, y_total_blocks
        """
        return self.x_total_blocks, self.y_total_blocks

    def set_block_size(self, block_width, block_height):

        self.block_padded_width = block_width
        self.block_padded_height = block_height

    def get_block_padded_size(self):
        """
        Gets size of bock with padding
        :return:
        """
        return self.block_padded_width, self.block_padded_height

    def set_block_bounds(self, block_top_left, block_bottom_right):
        """
        Internal use only, sets current block position and extent
        :param block_top_left:
        :param block_bottom_right:
        :return:
        """
        self.block_top_left = block_top_left
        self.block_bottom_right = block_bottom_right

    def get_block_coord_arrays(self):
        """
        Returns a numpy matrix of coordinates which represent real-world coordinates
        of each pixel on current coordinate arrays
        :return: x_block - x-coordinates of middle of each pixel, y_block - Y-coordinates of middle of each pixel
        """
        (top_left, bottom_right) = (self.block_top_left, self.block_bottom_right)
        (n_cols, n_rows) = self.get_block_padded_size()
        n_cols += 2 * self.padding
        n_rows += 2 * self.padding
        (x_res, y_res) = self.get_pixel_size()
        (row_ndx, col_ndx) = numpy.mgrid[0:n_rows, 0:n_cols]
        x_block = top_left.x - self.padding * x_res + x_res / 2.0 + col_ndx * x_res
        y_block = top_left.y + self.padding * y_res - y_res / 2.0 - row_ndx * y_res
        return x_block, y_block

    def set_current_block(self, x_block, y_block):
        """
        Internal use only , sets current padded block coordinates
        :param x_block:
        :param y_block:
        :return:
        """
        self.x_block = x_block
        self.y_block = y_block

    def get_block_count(self):
        """
        get total number of blocks
        :return:
        """
        return self.x_block, self.y_block

    def get_pixel_size(self):
        """
        Get resolution of current reference grid
        :return:
        """
        return self.reference_grid.x_res, self.reference_grid.y_res

    def get_pixel_row_column_block(self, x, y):
        """
        Gets row and column of pixel of current block, which correspond to x,y world coordinates
        :param x: world x coordinate
        :param y: world y coordinate
        :return: int block_row, int block_col or None, None of coordinate is out of bounds
        """
        transform = self.reference_grid.get_geotransform()
        img_row_col = utils.world_to_pixel(transform, x, y)

        # here they are, pixel coordinates of geographic point
        img_row = img_row_col.y
        img_col = img_row_col.x

        # now convert them to block coordinates
        block_start_row = self.y_block * self.block_size_y - self.padding
        block_start_col = self.x_block * self.block_size_x - self.padding

        # here they are
        block_row = int(img_row - block_start_row)
        block_col = int(img_col - block_start_col)

        # return None, if they are out of bounds
        if ((block_row < 0 or block_row > (self.block_size_y + 2 * self.padding)) or
                (block_col < 0 or block_col > (self.block_size_x + 2 * self.padding))):
            block_row = None
            block_col = None

        return block_row, block_col

    def get_pixel_column_row(self, x, y):
        """
        get pixel column and row for x (column) and y (row) pixel of current block
        :param x:
        :param y:
        :return:
        """
        col = self.x_block * self.block_size_x + x
        row = self.y_block * self.block_size_y + y
        return col, row

    def is_first_block(self):
        """
        :return: True if current block is first one
        """
        return self.x_block == 0 and self.y_block == 0

    def is_last_block(self):
        """
        :return: true, if current block is last one
        """
        return self.x_block == self.x_total_blocks - 1 and self.y_block == self.y_total_blocks - 1

    def get_filename_for(self, block):
        """
        Get dataset filename for specific block
        compares by block.__id__
        :param block: numpy array
        :return: string
        """
        (ds, fname) = self.block_lookup[id(block)]
        return fname

    def get_gdal_dataset_for(self, block):
        """
        Get gdal dataset for specific block
        block.__id__
        :param block: numpy array
        :return: gdal dataset
        """
        (ds, fname) = self.block_lookup[id(block)]
        return ds

    def get_gdal_band_for(self, block, band_index):
        """
        Get band at band_index of gdal dataset corresponding to current block
        :param block: numpy array
        :param band_index: integer (starting at 1)
        :return: ds.GetRasterBand(band_index)
        """
        ds = self.get_gdal_dataset_for(block)
        return ds.GetRasterBand(band_index)

    def get_no_data_value_for(self, block, band_index=1):
        """
        Get no data value for specific block and band
        :param block: numpy array
        :param band_index: integer, band index (starting at 1) (default - 1)
        :return:
        """
        ds = self.get_gdal_dataset_for(block)
        band_index = ds.GetRasterBand(band_index)
        novalue = band_index.GetNoDataValue()

        if novalue is not None:
            numpytype = utils.gdal_type_to_numpy_type(band_index.DataType)
            novalue = numpy.cast[numpytype](novalue)

        return novalue
