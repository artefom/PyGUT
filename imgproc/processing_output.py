import os
import math

from osgeo import gdal
import numpy

from . import errors
from . import constants
from . import utils


class ImageWriter(object):
    """
    Image writer class that supports writing blocks of data to output image
    designed to be used alongside with ImageReader
    """
    def __init__(self, filename, creationoptions=None, drivername=constants.DEFAULT_DRIVER_NAME, gdaldatatype=None,
                 nbands=None, info=None, firstblock=None, x_size=None, y_size=None, projection=None, transform=None,
                 block_size_x=None, block_size_y=None, padding=None):
        """
        :param filename: str, output filename
        :param creationoptions: creation options for gdal.Create list of strings, can't contain BLOCKXSIZE or BLOCKYSIZE
        :param drivername: str, name of gdal driver to write images with
        :param gdaldatatype: enum, gdal data type of output image
        :param nbands: int, number of bands
        :param info: info as passed from ImageReader
        :param firstblock: numpy array from which to determine number of bands and data type
        :param x_size: int, output image x size
        :param y_size: int, output image y size
        :param projection: projection WKT string
        :param transform: geo transformation list
        :param block_size_x: int, block x size
        :param block_size_y: int, block y size
        :param padding: padding, as passed from ImageReader
        """
        self.filename = filename
        noninfoitems = [x_size, y_size, transform, projection, block_size_x, block_size_y, padding]
        if info is None:
            # user did not pass info object, so make sure that he did pass all other objects
            if not all(not v is None for v in noninfoitems):
                msg = "make sure that either info object is not None or all ohters"
                raise errors.InvalidArgumentsError(msg)

            # get necessary data from parameters
            self.padding = padding
            self.block_size_x = block_size_x
            self.block_size_y = block_size_y
            self.x_total_blocks = int(math.ceil(float(x_size) / block_size_x))
            self.y_total_blocks = int(math.ceil(float(y_size) / block_size_y))

        else:
            # user actually did pass info object, retrieve information from it
            if any(not v is None for v in noninfoitems):
                msg = "If info object is passed, than all other parameters should be None"
                raise errors.InvalidArgumentsError(msg)

            # get necessary data from it
            (x_size, y_size) = info.get_total_size()
            transform = info.get_geotransform()
            projection = info.get_projection()
            (self.block_size_x, self.block_size_y) = info.get_block_size()
            self.padding = info.get_padding_size()
            (self.x_total_blocks, self.y_total_blocks) = info.get_total_blocks()

        if firstblock is None and not all(not v is None for v in [nbands, gdaldatatype]):
            msg = 'need first block or (nbands and gdal type) to be not None'
            raise errors.InvalidArgumentsError(msg)

        elif firstblock is not None and any(not v is None for v in [nbands, gdaldatatype]):
            msg = 'need first block or (nbands and gdal type) to be not None'
            raise errors.InvalidArgumentsError(msg)

        if firstblock is not None:
            # user passed first block, get number of bands and data type from it
            if len(firstblock.shape) != 3:
                raise errors.InvalidArrayShape(
                    "firshblock must have 3 dimensions, but have {}".format(repr(firstblock.shape)))

            (nbands, y, x) = firstblock.shape

            # map numpy type to gdal data type
            gdaldatatype = utils.numpy_type_to_gdal_type(firstblock.dtype)

        # make sure that creation options are not None
        if creationoptions is None:
            if drivername in constants.default_driver_parameters:
                creationoptions = constants.default_driver_parameters[drivername]
            else:
                creationoptions = []

        # add block size to creation options
        # at this moment this.window_size_x and y must be set
        creationoptions = self.make_creation_options(drivername, creationoptions)

        # self-explanatory
        self.delete_if_existing(filename)

        driver = gdal.GetDriverByName(drivername)
        self.ds = driver.Create(str(filename), x_size, y_size, nbands, gdaldatatype, creationoptions)
        if self.ds is None:
            msg = 'Failed to write result {} file'.format(filename)
            raise errors.OpenImageError(msg)

        # set projection and geotransform of newly created image
        self.ds.SetProjection(projection)
        self.ds.SetGeoTransform(transform)

        # set current block to first
        self.blockid = 0

        if firstblock is not None:
            self.write(firstblock)

    @staticmethod
    def delete_if_existing(filename):

        if os.path.exists(filename):

            using_exceptions = gdal.GetUseExceptions()
            if not using_exceptions:
                gdal.UseExceptions()

            try:
                ds = gdal.Open(str(filename))
            except RuntimeError:
                ds = None

            if ds is not None:

                drvr = ds.GetDriver()
                del ds

                # some datasets may contain additional files, so make gdal handle this
                drvr.Delete(filename)
            else:

                os.remove(filename)

            # and revive previous exceptions state
            if not using_exceptions:
                gdal.DontUseExceptions()

    def get_gdal_dataset(self):
        """
        :return: gdal dataset, current output dataset
        """
        return self.ds

    def get_current_block(self):
        """
        :return: int, index of current write block
        """
        return self.blockid

    def set_thematic(self):
        """
        set thematic layers
        :return: None
        """
        for i in range(1, self.ds.RasterCount + 1):
            band = self.ds.GetRasterBand(i)
            band.SetMetadataItem('LAYER_TYPE', 'thematic')

    def set_band_names(self, names):
        """
        set band names
        :param names: list of strings
        :return: None
        """
        bandindex = 1
        for name in names:
            bh = self.ds.GetRasterBand(bandindex)
            bh.Set_description(name)
            bandindex += 1

    def write(self, block):
        """
        write block to current block position
        :param block: numpy array, block to be written (with padding if any). usually get this one from user function
        :return:
        """
        y_block = self.blockid // self.x_total_blocks
        x_block = self.blockid % self.x_total_blocks

        xcoord = x_block * self.block_size_x
        ycoord = y_block * self.block_size_y

        self.write_at(block, xcoord, ycoord)

        self.blockid += 1

    def write_at(self, block, xcoord, ycoord):
        """
        Handle padding slicing and writing dataset
        :param block: numpy array
        :param xcoord: integer, top coordinate
        :param ycoord: integer, left coordinate
        :return:
        """

        # get block coordinates bottom right coordinates of block
        bottom_right_x_coord = xcoord + block.shape[-1] - self.padding * 2
        bottom_right_y_coord = ycoord + block.shape[-2] - self.padding * 2

        # block must match perfectly!
        if bottom_right_x_coord > self.ds.RasterXSize or bottom_right_y_coord > self.ds.RasterYSize:
            raise errors.ImageSizeError()

        # block must contain only 3 dimensions!
        if block.ndim != 3:
            raise errors.InvalidArgumentsError("block must contain only 3 dimensions!")

        # write to each band seperately
        for band in range(self.ds.RasterCount):
            bh = self.ds.GetRasterBand(band + 1)

            # slice off padding excess
            slice_bottom_most = block.shape[-2] - self.padding
            slice_right_most = block.shape[-1] - self.padding

            outblock = block[band, self.padding:slice_bottom_most, self.padding:slice_right_most]

            bh.WriteArray(outblock, xcoord, ycoord)

    def reset(self):

        self.blockid = 0

    def close(self):

        self.ds.FlushCache()
        del self.ds
        self.ds = None

    def make_creation_options(self, drivername, creationoptions):
        """
        Mekeup creation options for file writing
        More specially, add block sizes if drivername is 'GTiff'
        :param drivername: string, name of gdal driver
        :param creationoptions: list of strings - options for gdal.Create()
        :return:
        """
        new_creationoptions = creationoptions

        if drivername == 'GTiff':

            img_block_x = self.block_size_x
            img_block_y = self.block_size_y

            tiff_block_x = img_block_x
            tiff_block_y = img_block_y

            def is_power_of2(n):
                return ((n - 1) & n) == 0

            if not (is_power_of2(tiff_block_x) and is_power_of2(tiff_block_y)):
                msg = "GTiff block size must be power of 2. Have: {} ".format((tiff_block_x, tiff_block_y))
                raise errors.OpenImageError(msg)

            new_creationoptions.append('BLOCKXSIZE={}'.format(self.block_size_x))
            new_creationoptions.append('BLOCKYSIZE={}'.format(self.block_size_y))

        return new_creationoptions
