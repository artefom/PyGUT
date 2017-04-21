# Artyom Fomenko BSE142

import sys

from osgeo import ogr
from osgeo import gdal
import numpy

from . import processing_output
from . import processing_input
from . import errors
from . import constants


class FilenameContainer(object):
    """
    Class containing information about input objects
    Filenames are values of class fields as string or list of strings
    (see __dict__)
    """
    def __len__(self):
        return len(self.__dict__.keys())


class BlockContainer(object):
    """
    Class containing information about processing blocks.
    Blocks or list of blocks are values of class fields
    see __dict__
    """
    def __len__(self):
        return len(self.__dict__.keys())


class ArgumentsContainer(object):
    """
    Class containing arguments passed to user function.
    Example::
    c = ArgumentsContainer()
    c.arg1 = 10
    c.arg2 = "hi"

    apply(inputs,outputs,c)
    """
    pass


class CalculationParams(object):
    """
    Class, representing parameters for applier function.
    Parameters :
        self.logstream - stream to write log into
        self.drivername - gdal driver name, used for file writing (default: GTiff)
        self.padding - number of padding pixels of block
        self.window_size_x - x size of window (use power of 2)
        self.window_size_y - y size of windows (use power of 2)
        self.footprint - result trimmed to include union, intersection, or custom footpring
        self.reference_image - image, from what derive result rastergrid
            (used with self.footprint set to REFERENCE)
        self.reference_rastgrid - raster_grid, that determines result area
            (used with self.footprint set to REFERENCE)
        self.creationoptions - parameters for gdal creation script
        self.band_names - list of layer names for result
        self.tempdir - temporary directory for resampled images
        self.resample_method - default method for image resampling
    """

    def __init__(self):
        self.logstream = sys.stdout  # stream to write log into
        self.drivername = constants.DEFAULT_DRIVER_NAME  # gdal driver name, used for file writing (default: GTiff)
        self.padding = constants.DEFAULT_PADDING  # number of padding pixels of block
        self.window_size_x = constants.DEFAULT_WINDOW_SIZE_X  # x size of window (use power of 2)
        self.window_size_y = constants.DEFAULT_WINDOW_SIZE_Y  # y size of windows (use power of 2)
        self.footprint = constants.DEFAULT_FOOTPRINT  # Result trimmed to include union, intersection, or custom footpring
        self.reference_image = None  # image, from what derive result rastergrid
        self.reference_rastgrid = None  # raster_grid, that determines result area
        self.creationoptions = None  # parameters for gdal creation script
        self.thematic = None  # determines weather the result should be thematic
        self.band_names = None  # list of layer names for result
        self.layer_selection = None  # list of selected layers
        self.tempdir = '.'  # temporary directory for resampled images
        self.resample_method = constants.DEFAULT_RESAMPLEMETHOD

        self.options_by_image = {}

    def set_option_for_imagename(self, option, imagename, value):
        """
        Set unique options for image.
        These options are used on image read and write.
        for read, layer_selection and resamplemethod is derived from options
        for write, such things as drivername, thematic, band names, etc..

        Duplicates fields of current class, priority will be given to options specific to image

        :param option: option to be set
        :param imagename: name of image for current option (must me in inputs our outputs containers)
                           None to set for all images
        :param value: value of option
        :return: None
        """
        if imagename is None:
            setattr(self, option, value)
        else:
            if not option in self.options_by_image:
                self.options_by_image[option] = {}
            self.options_by_image[option][imagename] = value

    def get_option_for_imagename(self, option, imagename):
        """
        Get option specific for image. see set_option_for_imagename()
        :param option:
        :param imagename:
        :return:
        """
        value = getattr(self, option)
        if option in self.options_by_image:
            if imagename in self.options_by_image[option]:
                value = self.options_by_image[option][imagename]
        return value

    def set_logstream(self, logstream):
        """
        Set stream which to write debuf info into.
        :param logstream: file-like object (sys.stdout for default)
        :return: None
        """
        self.logstream = logstream

    def set_padding(self, padding):
        """
        When splitting image into blocks, they may overlap by some pixels.
        padding is number of pixels which overlap for each block
        :param padding: integer
        :return: None
        """
        self.padding = padding

    def set_result_driver_name(self, drivername, imagename=None):
        """
        Set GDAL driver to write result. GTiff by default
        :param drivername: One of gdal drivers (see http://www.gdal.org/formats_list.html)
        :param imagename: string name of output image or None to set for all
        :return: None
        """
        self.set_option_for_imagename('drivername', imagename, drivername)

    def set_window_x_size(self, window_size_x):
        """
        Set x size of processing block
        :param window_size_x: int
        :return: None
        """
        self.window_size_x = window_size_x

    def set_window_y_size(self, window_size_y):
        """
        Set y size of processing block
        :param window_size_y: int
        :return: None
        """
        self.window_size_y = window_size_y

    def set_footprint_type(self, footprint):
        """
        Type of resulting footprint (area of result)
        It may be either FP_INTERSECTION for intersection of input areas
        or FP_UNION for union of input areas
        or FP_REFERENCE to derive output footprint from reference image
        (see set_reference)

        :param footprint: enum (constants.FP_INTERSECTION by default)
        :return: None
        """
        self.footprint = footprint

    def set_creation_options(self, creationoptions, imagename=None):
        """
        Set some arbitrary creation options for image.
        :param creationoptions: list of strings
        :param imagename: name of image
        :return:
        """
        self.set_option_for_imagename('creationoptions', imagename, creationoptions)

    def set_thematic(self, thematic_flag, imagename=None):
        """
        Weather output should be thematic or not
        :param thematic_flag:
        :param imagename:
        :return:
        """
        self.set_option_for_imagename('thematic', imagename, thematic_flag)

    def set_band_names(self, band_names, imagename=None):
        """
        Set band names for specific image, or all images
        :param band_names: list of strings?
        :param imagename:
        :return:
        """
        self.set_option_for_imagename('band_names', imagename, band_names)

    def set_tempdir(self, tempdir):
        """
        Set temporary dir where to store all temporary files when resampling
        All temporary files are cleaned shortly after they used
        :param tempdir: relative path to directory. By default - "."
        :return: None
        """
        self.tempdir = tempdir

    def set_resample_method(self, resample_method, imagename=None):
        """
        Set resample method for input images
        :param resample_method: one of following: 'bilinear', 'near', 'cubicspline',
        'cubic', 'lanczos'. see gdalwarp for more info
        :param imagename: specific image or None for all images
        :return: None
        """
        self.set_option_for_imagename('resample_method', imagename, resample_method)

    def workout_resample_dictionary(self, image_dict):
        """
        Internal function. used to make resample dictionary for each image
        :param image_dict: dictionary of input images wher keys are their names and values are paths
        :return: dict, containing where keys are image names and values are resample methods
        """
        d = {}
        imagenamelist = image_dict.keys()
        for name in imagenamelist:
            method = self.get_option_for_imagename('resample_method', name)
            if isinstance(image_dict[name], list):
                d[name] = [method] * len(image_dict[name])
            else:
                d[name] = method
        return d

    def select_input_image_layers(self, layer_selection, imagename=None):

        self.set_option_for_imagename('layer_selection', imagename, layer_selection)

def calculate(custom_function, input_files, output_files, custom_fun_params=None, parameters=None):
    """
    Applies 'custom_function' to given input and output files

    usage example:

    input_files = FilenameContainer()
    input_files.img1 = 'image1.tiff'
    input_files.img2 = 'image2.tiff'
    output_files = FilenameContainer()
    output_files.output_image = 'output.tiff'

    :param custom_function: function with signature foo(inputs,outputs,info)->None
    :param input_files: FilenameContainer, containing image names as values of attributes
    :param output_files: FilenameContainer, containing image names as values of attributes
    :param custom_fun_params: additional parameters passed to user function
    :param parameters: CalculationParams which control calculationProcess
    :return:
    """
    if parameters is None:
        parameters = CalculationParams()

    imagefiles = images_validate(input_files)
    input_image_layer_selection = get_input_layer_selection(imagefiles, parameters)
    reader = processing_input.ImageReader(imagefiles.__dict__,
                                          parameters.footprint, parameters.window_size_x, parameters.window_size_y,
                                          parameters.padding, logstream=parameters.logstream,
                                          layer_selection=input_image_layer_selection)

    resample_inputs(imagefiles, parameters, reader)

    # dictionary of writers for each image
    writerdict = {}

    done = False
    iterator = reader.__iter__()
    while not done:

        try:
            info, blockdict = iterator.__next__()
        except StopIteration:
            done = True
            break

        input_blocks = BlockContainer()
        input_blocks.__dict__.update(blockdict)

        params = get_function_params(info, input_blocks, custom_fun_params)
        custom_function(*params)

        result_blocks = params[1]
        write_result_blocks(writerdict, output_files, result_blocks,
                            parameters, info)

    close_output_images(writerdict, output_files, parameters)


def close_output_images(writerdict, output_files, parameters):
    """
    internal function.
    called by calculate() to close all resulting image files
    :param writerdict: dictionary of Imagewriters
    :param output_files: FilenameContainer, names of output files
    :param parameters: CalculationParams
    :return:
    """
    for name in output_files.__dict__.keys():
        writer = writerdict[name]
        if isinstance(writer, list):
            writer_list = writer
        else:
            writer_list = [writer]

        for single_writer in writer_list:
            single_writer.close()


def resample_inputs(input_files, parameters, reader):
    """
    Handles input resampling based on parameters.

    In come way, this function is a way to control reader's behaviour
    through parameters argument

    calls reader.resample_inputs with resample methods derived from parameters
    based on parameters.reference_image or parameters_rastgrid
    :param input_files: FilenameContainer of input files
    :param parameters: CalculationParams
    :param reader: ImageReader
    :return:
    """
    if parameters.reference_image is not None:
        resample_dict = parameters.workout_resample_dictionary(input_files.__dict__)
        reader.resample_inputs(refpath=parameters.reference_image, tempdir=parameters.tempdir,
                               resampling_method=resample_dict)
    elif parameters.reference_rastgrid is not None:
        resample_dict = parameters.workout_resample_dictionary(input_files.__dict__)
        reader.resample_inputs(ref_rastgrid=parameters.reference_rastgrid,
                               tempdir=parameters.tempdir,
                               resampling_method=resample_dict)


def write_result_blocks(writerdict, output_files, result_blocks, parameters, info):
    """
    Writes resulting blocks, created by custom user function.
    Also, creates writers on first call.
    :param writerdict: dictionary of writers for each file
    :param output_files: FilenameContainer of output files
    :param result_blocks: BlockContainer of created blocks
    :param parameters: CalcilationParams
    :param info:
    :return:
    """
    for name in output_files.__dict__.keys():
        if name not in result_blocks.__dict__:
            msg = 'Result key {} not found in result blocks'.format(name)
            raise errors.InvalidKey(msg)

        outblock = result_blocks.__dict__[name]
        outfile_name = getattr(output_files, name)
        if name not in writerdict:
            # result writers are not created yet, so create them
            if isinstance(outfile_name, list):
                # create list of writers for list of ouput files
                writerdict[name] = []
                num_files = len(outfile_name)
                if len(outblock) != num_files:
                    raise errors.ListLengthError(("Result '{}' writes {} files, " +
                                                  "but only {} blocks given").format(name, num_files, len(outblock)))
                for i in range(num_files):
                    filename = outfile_name[i]
                    writer = processing_output.ImageWriter(filename,
                                                           creationoptions=parameters.get_option_for_imagename(
                                                               'creationoptions', name),
                                                           drivername=parameters.get_option_for_imagename('drivername',
                                                                                                          name),
                                                           info=info, firstblock=outblock[i])
                    writerdict[name].append(writer)
                    if parameters.get_option_for_imagename('thematic', name):
                        writer.set_thematic()

                    band_names = parameters.get_option_for_imagename('band_names', name)
                    if band_names is not None:
                        writer.set_band_names(band_names)
            else:
                # or create single writer for single output filename
                writer = processing_output.ImageWriter(outfile_name,
                                                       creationoptions=parameters.get_option_for_imagename(
                                                           'creationoptions',
                                                           name),
                                                       drivername=parameters.get_option_for_imagename('drivername',
                                                                                                      name), info=info,
                                                       firstblock=outblock)
                writerdict[name] = writer
                if parameters.get_option_for_imagename('thematic', name):
                    writer.set_thematic()

                band_names = parameters.get_option_for_imagename('band_names', name)
                if band_names is not None:
                    writer.set_band_names(band_names)
        else:
            # result writes have already been created, select correct one

            # if output file is list of strings:
            if isinstance(outfile_name, list):
                num_files = len(outfile_name)
                if len(outblock) != num_files:
                    raise errors.ListLengthError(("Result '{}' writes {} files, " +
                                                  "but only {} blocks given").format(name, num_files, len(outblock)))
                for i in range(num_files):
                    # select corresponding writer and write block
                    writerdict[name][i].write(outblock[i])
            else:
            #outptut file is just a string, so write it strait forward
                writerdict[name].write(outblock)


def images_validate(infiles):
    """
    takes FilenameContainer and ensures that it's fields are strings
    or lists of strings only, tries to open them (with try_open)
    and yields if gdal had thrown any exceptions.
    :param infiles: FilenameContainer
    :return: newly created filename container with all valid files
    """
    imagefiles = FilenameContainer()

    name_list = sorted(infiles.__dict__.keys())
    for name in name_list:
        filename = getattr(infiles, name)
        if isinstance(filename, str):
            test_filename = filename
        elif isinstance(filename, list):
            # got list of filenames, just test first one
            test_filename = filename[0]
        else:
            # filename is some other unknown shit, skip it.
            # or maybe it is not even a filename
            test_filename = None
        if try_open(test_filename):
            setattr(imagefiles, name, filename)
        else:
            raise errors.OpenFileError("Failed to open file '{}'".format(test_filename))

    return imagefiles


def try_open(filename):
    """
    tries to open image file with gdal.Open and throws an exception
    if something went wrong.
    :param filename:
    :return:
    """
    # make sure gdal yields all exceptions
    using_exceptions = False
    if hasattr(gdal, 'GetUseExceptions'):
        using_exceptions = gdal.GetUseExceptions()
    gdal.UseExceptions()
    try:
        ds = gdal.Open(filename)
    except Exception:
        ds = None
    opens_o_k = (ds is not None)

    # and revive previous state
    if not using_exceptions:
        gdal.DontUseExceptions()
    return opens_o_k


def get_input_layer_selection(imagefiles, controls):
    """
    Get selection of input layers for imagefiles
    :param imagefiles: FilenameContainer
    :param controls: CalculationParams
    :return:
    """
    layer_selection = dict()
    for name in imagefiles.__dict__.keys():
        layer_selection[name] = controls.get_option_for_imagename('layer_selection', name)
    return layer_selection

def get_function_params(info, inputs, otherargs=None):
    """
    Get params for user custom function
    :param info: current calculation process info
    :param inputs: input filenames
    :param otherargs: custom user arguments
    :return: tuple (inputs,results,info,*otherargs) of parameters for user function
    """
    if info.logstream is None:
        info.logstream = constants.DEFAULT_LOGSTREAM
    results = BlockContainer()
    params = (inputs, results, info)
    if otherargs is not None:
        params += (otherargs,)

    return params