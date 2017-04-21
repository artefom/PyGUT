import sys
import copy
import numpy
from osgeo import gdal
from . import utils
from . import inputcontainer
from . import readerinfo
from . import errors
from . import constants


class ImageIterator(object):
    """
    Iterator of Imagereader class
    Do not create this class
    """
    def __init__(self, reader):
        self.reader = reader
        self.nblock = 0

    def __iter__(self):
        return self

    def __next__(self):
        # catch out of bounds exception of read_block and re-raise StopIteration exception
        # supporting python standards
        try:
            # get tuple of curent block
            return_tuple = self.reader.read_block(self.nblock)
        except errors.ImageSizeError:
            raise StopIteration()

        # look at next block on next iteration
        self.nblock += 1

        return return_tuple


class ImageReader(object):
    """
    Class supporting reading input files block by block and responsible for input resampling
    """
    def __init__(self, image_container,
                 footprint=constants.DEFAULT_FOOTPRINT,
                 window_size_x=constants.DEFAULT_WINDOW_SIZE_X,
                 window_size_y=constants.DEFAULT_WINDOW_SIZE_Y,
                 padding=constants.DEFAULT_PADDING,
                 logstream=sys.stdout,
                 layer_selection=None):
        """
        Constructor
        :param image_container: FilenameContainer, containing input images
        :param footprint: enum default: FP_INTERSECTION
            working area is derived either as intersection of input images
            (FP_INTERSECTION)
            or union
            (FP_UNION)
            or based on refrence image
            (FP_REFERENCE)
        :param window_size_x:
        :param window_size_y:
        :param padding:
        :param logstream:
        :param layer_selection:
        """
        self.image_container = image_container

        if isinstance(image_container, dict):
            # passed dictionary of files
            image_list = []
            self.layer_selection_list = []
            for name in image_container.keys():
                filename = image_container[name]
                if isinstance(filename, str):

                    image_list.append(filename)
                else:
                    msg = "image container must contain strings. Got {}".format(type(filename))
                    raise errors.InvalidArgumentsError(msg)

                # get current selection for layer
                if layer_selection is None:
                    this_layer_selection = None
                else:
                    this_layer_selection = layer_selection[name]

                self.layer_selection_list.append(this_layer_selection)


        elif isinstance(image_container, str):
            # passed single file
            image_list = [image_container]
            if layer_selection is not None:
                self.layer_selection_list = [layer_selection]
            else:
                self.layer_selection_list = [None]
        else:
            # passed list of images
            image_list = image_container
            if layer_selection is not None:
                self.layer_selection_list = layer_selection
            else:
                self.layer_selection_list = [None for _ in image_list]

        # create container for input images
        self.inputs = inputcontainer.InputContainer(image_list, logstream=logstream)

        # remember passed information
        self.footprint = footprint
        self.window_size_x = window_size_x
        self.window_size_y = window_size_y
        self.padding = padding
        self.logstream = logstream

        # these are set later during prepare method
        self.working_grid = None
        self.info = None

    def __len__(self):
        """
        :return: total number of blocks
        """
        if self.info is None:
            self.prepare()

        (x_total_blocks, y_total_blocks) = self.info.get_total_blocks()

        return x_total_blocks * y_total_blocks

    def __getitem__(self, key):
        """
        Get block by key index
        :param key: integer, block index
        :return: ReaderInfo, BlockContainer
        """
        if self.info is None:
            self.prepare()

        if key < 0:
            # if key is negative, offset it from end
            (x_total_blocks, y_total_blocks) = self.info.get_total_blocks()

            key = (x_total_blocks * y_total_blocks) + key
            if key < 0: # key is still negative, not enough blocks
                raise KeyError()

        try:
            # try reading block y key
            return_tuple = self.read_block(key)
        except errors.ImageSizeError:
            # oops, wrong index anyway
            raise KeyError()

        # return this crap info, BockContainer
        return return_tuple

    def __iter__(self):
        """
        Get iterator for current inputs
        :return:
        """
        if self.info is None:
            self.prepare()

        return ImageIterator(self)

    def resample_inputs(self, resampling_method="near", refpath=None, refgeotrans=None,
                        refproj=None, ref_n_cols=None, ref_n_rows=None, ref_rastgrid=None,
                        tempdir='.'):

        """
        Refactor all inputs to reference grid
        :param resampling_method: 'near' by default see gdalwarp for other options
        :param refpath: path to reference file ( from which extract resulting area )
        :param refgeotrans: reference geoposition. must be passed if refpath is None and ref_rastgrid is None
        :param refproj: refence projection. must be passed if refpath is None and ref_rastgrid is None
        :param ref_n_cols: number of culumns. must be passed if refpath is None and ref_rastergrid is None
        :param ref_n_rows: number of rows. must be passed it refpath is None and ref_rastergrid is None
        :param ref_rastgrid: raster grid to be used as reference. must be passed if refpath is None
        :param tempdir:
        :return:
        """
        self.inputs.set_reference(refpath, refproj, refgeotrans, ref_n_cols, ref_n_rows, ref_rastgrid)

        if isinstance(resampling_method, str):
            # resample method passed as single string, construct list of resample methods for each input by copying
            resample_methods_lst = [resampling_method] * len(self.inputs)
        elif isinstance(resampling_method, dict):
            # resample methods passed as dictionary
            if not isinstance(self.image_container, dict):
                msg = 'can only pass resamplemethod as dict if dict was initially passed to constructor'
                raise errors.InvalidArgumentsError(msg)
            elif sorted(self.image_container.keys()) != sorted(resampling_method.keys()):
                #
                msg = ('dictionary keys in input image_container must match those in dictionary passed to constructor' +
                       'image keys = {}, resamplemethod keys = {}').format(self.image_container.keys(),
                                                                           resampling_method.keys())
                raise errors.InvalidArgumentsError(msg)
            else:
                # resamplemethod passed as dictionary, everything is OK, all checks passed

                # construct resample_methods_lst from dictionary
                resample_methods_lst = []

                # iterate through keys and add values to resample_methods_lst
                for name in resampling_method.keys():
                    method = resampling_method[name]
                    if isinstance(method, list):
                        # if list encountered, extend resample_methods_lst by it
                        resample_methods_lst.extend(method)
                    elif isinstance(method, str):
                        # if single string, add it ti resample_methods_lst
                        resample_methods_lst.append(method)
                    else:
                        msg = "resample_methods_lst can contain only str or list. Have: {}".format(type(method))
                        raise errors.InvalidArgumentsError(msg)

        else:
            # passed list of resampling methods
            # check if it size is valid
            if len(resampling_method) != len(self.inputs):
                msg = 'resampling_methods length must equal to list passed to constructor'
                raise errors.InvalidArgumentsError(msg)
            resample_methods_lst = resampling_method

        try:
            # all resample methods are figured out, make resampling.
            self.inputs.resample_all_to_reference(self.footprint, resample_methods_lst, tempdir)
        finally:
            # even if everything is failed, cleanup all temporary files
            self.inputs.cleanup()

    def prepare(self, working_grid=None):
        """
        Preparing to read from input images,
        :param working_grid:
        :return:
        """

        # check that all images match
        if not self.inputs.is_all_match():
            msg = "Images don't match, please, resample them"
            raise errors.ResampleRequiedError(msg)

        # set reference grid to supplied or retrieve it from footprint type
        self.working_grid = self.inputs.find_working_region(self.footprint) if working_grid is None else working_grid

        # create reader info
        # later, a shallow copy of this info will be pased to user function
        self.info = readerinfo.ReaderInfo(self.working_grid,
                                          self.window_size_x, self.window_size_y, self.padding, self.logstream)

    def get_block_position(self,nblock):
        """
        convert block index into it's x and y position
        :param nblock: integer, block index
        :return: x: int, y: int
        """
        (x_total_blocks, y_total_blocks) = self.info.get_total_blocks()

        # raise error if block index is out of bounds
        if nblock >= (x_total_blocks * y_total_blocks):
            raise errors.ImageSizeError()

        y_block = nblock // x_total_blocks
        x_block = nblock % x_total_blocks

        return x_block, y_block

    def read_block(self, nblock):
        """
        Return specific block by it's index
        nblock will be converted into block's x and y
        Called when iterating over current class instance
        :param nblock: block index
        :return: ReaderInfo, BlockContainer
        """

        # again, make sure everything is prepared
        if self.info is None:
            self.prepare()

        # make shallow copy of info to pass ot to user
        info = copy.copy(self.info)

        # convert nblock to x_block and y_block
        x_size, y_size = info.get_total_size()

        x_block, y_block = self.get_block_position(nblock)

        # update current block position at info
        info.set_current_block(x_block, y_block)

        # get coordinates of top left pixel of current block
        xcoord = x_block * self.window_size_x
        ycoord = y_block * self.window_size_y

        # get world coordinates of top left of top-left pixel of current block
        block_top_left = utils.pixel_to_world(info.get_geotransform(), xcoord, ycoord)

        # get coordinates of bottom+1 left+1 pixel of current block
        n_block_bottom_x = min(((x_block + 1) * self.window_size_x),x_size)
        n_block_bottom_y = min(((y_block + 1) * self.window_size_y),y_size)

        # get world coordinates of bottom left of bottom-left pixel of block
        block_bottom_right = utils.pixel_to_world(info.get_geotransform(), n_block_bottom_x, n_block_bottom_y)

        # set block bounds
        info.set_block_bounds(block_top_left, block_bottom_right)

        # calculate block dimensions
        block_width = n_block_bottom_x - xcoord
        block_height = n_block_bottom_y - ycoord

        # set block size
        info.set_block_size(block_width, block_height)

        # begin creating tuple passed to user function
        block_list = []

        try:

            for layer_selection,(image, ds, rastgrid, null_val_list, datatype) in zip(self.layer_selection_list, self.inputs ):
                top_left = utils.world_to_pixel(rastgrid.get_geotransform(), block_top_left.x, block_top_left.y)

                # read block with specified overlap
                block = self.read_block_with_padding(ds, int(round(top_left.x)), int(round(top_left.y)), block_width, block_height,
                                                     datatype, padding=self.padding, null_val_list=null_val_list,
                                                     layer_selection=layer_selection)

                # append block to our tuple
                block_list.append(block)

                info.set_block_dataset(block, ds, image)

        finally:
            # cleanup temporary files
            self.inputs.cleanup()

        # if user passed dictionary, transform block list into dictionary
        if isinstance(self.image_container, dict):

            # block dicionary by itself
            block_dict = {}
            i = 0
            for name in self.image_container.keys():
                filename = self.image_container[name]
                if isinstance(filename, list):
                    list_len = len(filename)
                    # create list of blocks for list entry
                    block_dict[name] = []
                    # and fill it
                    for j in range(list_len):
                        block_dict[name].append(block_list[i])
                        i += 1
                elif isinstance(filename, str):
                    # assign block to corresponding image name key
                    block_dict[name] = block_list[i]
                    i += 1

            # reutrn value is block dictionary
            block_container = block_dict

        elif isinstance(self.image_container, str):
            # user passed single string, so
            # return value is first entry of block list
            block_container = block_list[0]

        else:
            # seems like user passed list, so return a list
            block_container = tuple(block_list)

        return info, block_container

    @staticmethod
    def read_block_with_padding(ds, xoff, yoff, x_size, y_size, datatype, padding=0, null_val_list=None,
                                layer_selection=None):
        """
        Function that reads block of data from ds with padding.
        Function returns valid numpy array with size x_size + 2*padding, y_size+2*padding in any case
        This function will return valid numpy array even if the area completely outsize ds image
        Values that lie outside image boundaries will be filled with values in null_val_list.

        :param ds: dataset from which to read array using ReadAsArray
        :param xoff: x offset of block
        :param yoff: y offset of block
        :param x_size: x size of block
        :param y_size: y size of block
        :param datatype: data type of returned numpy array
        :param padding: number of pixels to 'pad' around block
        :param null_val_list: list of values to fill no data pixels with.
            layer with index layer_selection[0] will be filled with null_val_list[0]
            layer with index layer_selection[1] will be filled with null_val_list[1]
            and so on...
        :param layer_selection: indecies of selected layers(bands) first band has index 1.
        :return:
        """

        # generate layer selection if it's none
        if layer_selection is None:
            layer_selection = [i + 1 for i in range(ds.RasterCount)]
        n_layers = len(layer_selection)

        x_size_padding = x_size + 2 * padding
        y_size_padding = y_size + 2 * padding
        out_block_shape = (n_layers, y_size_padding, x_size_padding)

        block_overlap = numpy.zeros(out_block_shape, dtype=datatype)
        if null_val_list is not None and len(null_val_list) > 0:
            # if any of null values are None, assign them to 0
            # fill val list now contains valid fill values
            fill_val_list = [null_val for null_val in null_val_list]
            for i in range(len(fill_val_list)):
                if fill_val_list[i] is None:
                    fill_val_list[i] = 0

            # fill block with null values
            if len(out_block_shape) == 2:
                # fill whole block if number of it's dimensions is 2
                block_overlap.fill(fill_val_list[0])
            else:
                # or band by band with corresponding values
                for i in range(n_layers):
                    block_overlap[i].fill(fill_val_list[layer_selection[i] - 1])

        def clamp(x,v_min,v_max): return min(max(x,v_min),v_max)

        # new offset with padding
        xoff_padding = xoff - padding
        yoff_padding = yoff - padding

        # clamp these values
        xoff_padding_clamped = clamp(xoff_padding,0,ds.RasterXSize)
        x_right_padding_clamped = clamp( xoff_padding + x_size_padding, 0, ds.RasterXSize )
        x_size_padding_clamped = x_right_padding_clamped - xoff_padding_clamped

        # also, clamp y values
        yoff_padding_clamped = clamp(yoff_padding, 0, ds.RasterYSize)
        y_bottom_padding_clamped = clamp( yoff_padding + y_size_padding, 0, ds.RasterYSize )
        y_size_padding_clamped = y_bottom_padding_clamped - yoff_padding_clamped

        # estimate number of pixels that are not read from file
        # number of block left pixels that are not read from left border
        not_read_left = xoff_padding_clamped - xoff_padding
        # number of pixels that are not read from right border
        not_read_right = x_size_padding - not_read_left - x_size_padding_clamped
        # number of pixels that are not read from top
        not_read_top = yoff_padding_clamped - yoff_padding
        # number of pixels that are not read from bottom
        not_read_bottom = y_size_padding - not_read_top - y_size_padding_clamped

        # get slices of area to be read
        slice_right = x_size_padding - not_read_right
        slice_bottom = y_size_padding - not_read_bottom

        if x_size_padding_clamped > 0 and y_size_padding_clamped > 0:

            # generate slice object
            image_slice = (slice(not_read_top, slice_bottom), slice(not_read_left, slice_right))

            # assign read array to sliced area of block_overlap for each band
            for i in range(n_layers):
                band = ds.GetRasterBand(layer_selection[i])
                block_overlap[i][image_slice] = band.ReadAsArray(xoff_padding_clamped, yoff_padding_clamped,
                                                                   x_size_padding_clamped, y_size_padding_clamped)

        return block_overlap

    def close(self):

        self.inputs.close()
