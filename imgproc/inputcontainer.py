# file manipulations
import os

# logstream here
import sys

# library to execute gdal warp
import subprocess
from osgeo import gdal
from . import rastergrid
from . import errors
from . import utils
from . import constants

class InputIterator(object):
    """
    Allows iteration over InputContainer
    See python documentation on iterators for more

    :return
    image name (str),
    gdal dataset,
    raster grid,
    null value replacement,
    data type
    """

    def __next__(self):

        # Check if iteration is finished and raise StopIteration()
        # done acoording to python standards
        if self.index >= len(self.collection.dataset_list):
            raise StopIteration()

        # retrieve all necessary information based on current index
        image = self.collection.image_list[self.index]
        ds = self.collection.dataset_list[self.index]
        rastgrid = self.collection.rastgrid_list[self.index]
        nullvals = self.collection.null_val_list[self.index]
        datatype = self.collection.data_type_list[self.index]

        # increase current position by one
        self.index += 1

        # and return it
        return image, ds, rastgrid, nullvals, datatype


    def __init__(self, collection):
        # collection - instance of InputCollection
        self.collection = collection
        self.index = 0

    def __iter__(self):
        return self

class InputContainer(object):
    """
    Class designed to keep information about all inputs.

    Gets list of filenames, opens them and retrieves their values
    """
    def __init__(self, image_list, logstream=sys.stdout):

        self.logstream = logstream

        # initialize storage for values
        self.image_list = [] # list of image paths
        self.dataset_list = [] # datasets
        self.rastgrid_list = [] # list of raster grids derived from datasets
        self.null_val_list = [] # list of lists, null values for each band
        self.data_type_list = [] # list of numpy types
        self.temporary_files = [] # list of temporary files

        # for every image
        for image in image_list:

            # open it
            ds = gdal.Open(str(image))
            if ds is None:
                msg = 'Unable to Open {}'.format(image)
                raise errors.OpenImageError(msg)

            # create raster grid from dataset
            rast_grid = self.make_rast_grid_from_dataset(ds)

            # get null values for each band
            ds_null_val_list = []
            for i in range(ds.RasterCount):
                null_val = ds.GetRasterBand(i + 1).GetNoDataValue()
                ds_null_val_list.append(null_val)

            # get data type from first band (bands begin with 1)
            gdaldatatype = ds.GetRasterBand(1).DataType
            numpytype = utils.gdal_type_to_numpy_type(gdaldatatype)

            # store all retrieved values
            self.image_list.append(image)
            self.dataset_list.append(ds)
            self.rastgrid_list.append(rast_grid)
            self.null_val_list.append(ds_null_val_list)
            self.data_type_list.append(numpytype)

        # set reference datase as first file
        self.reference_rast_grid = self.rastgrid_list[0]

    def __del__(self):
        # on destruction, remove all temporary files
        self.cleanup()

    def cleanup(self):
        """
        Remove all temporary files
        :return: None
        """
        for f in self.temporary_files:
            if os.path.exists(f):
                os.remove(f)
        self.temporary_files = []

    def close(self):
        """
        close all opened datasets
        (causes removal of all temoprary files)
        :return: None
        """
        for ds in self.dataset_list:
            del ds
        self.dataset_list = []
        self.cleanup()

    def __getitem__(self, key):
        """
        Get stored parameters buy key index
        :param key: ineger, key
        :return: image, dataset, raster grid, null value, data type
        """
        image = self.image_list[key]
        ds = self.dataset_list[key]
        rastgrid = self.rastgrid_list[key]
        null_val_list = self.null_val_list[key]
        datatype = self.data_type_list[key]
        return image, ds, rastgrid, null_val_list, datatype

    def __len__(self):
        """
        get number of opened datasets
        :return: len(self.dataset_list)
        """
        return len(self.dataset_list)

    def __iter__(self):
        """
        Supports iteration!
        :return: InputIterator(self)
        """
        return InputIterator(self)

    def set_reference(self, file_path=None, projection=None, geo_transformation=None, ref_n_cols=None, ref_n_rows=None,
                      ref_rastgrid=None):
        """
        Retrieve reference grid to resample all inputs into it.

        You can provide either:
            file_path - path to dataset, from which raster grid will be recovered
        or
            projection, transformation, number of columns and rows
        or
            reference raster grid

        :param file_path: path to open as gdal dataset
        :param projection: projection wkt - string
        :param geo_transformation: transformation values - list
        :param ref_n_cols: number of columns - int
        :param ref_n_rows: number of rows - int
        :param ref_rastgrid: raster grid
        :return: None
        """

        # provied path to dataset!
        if file_path is not None:
            ds = gdal.Open(file_path)
            if ds is None:
                msg = 'Unable to Open {}'.format(file_path)
                raise errors.OpenImageError(msg)

            self.reference_rast_grid = self.make_rast_grid_from_dataset(ds)

            del ds

        # provided transformations!
        elif geo_transformation is not None and projection is not None and ref_n_cols is not None and ref_n_rows is not None:

            proj = projection
            rast_grid = rastergrid.RasterGrid(projection=proj, geotransform=geo_transformation,
                                              row_number=ref_n_rows, column_number=ref_n_cols)
            self.reference_rast_grid = rast_grid

        # provided raster grid
        elif ref_rastgrid is not None:

            ref_rastgrid.projection = ref_rastgrid.projection
            self.reference_rast_grid = ref_rastgrid
        else:
            msg = 'file_path or ref_rastgrid or all other parameters must be provided'
            raise errors.InvalidArgumentsError(msg)

    def resample_to_reference(self, ds, null_val_list, working_region, resamplemethod, tempdir='.'):
        """
        Resamples dataset to working_region
        :param ds: Dataset to be resampled
        :param null_val_list: value with which to fill empty regions
        :param working_region: raster grid defn
        :param resamplemethod: method to resample with
        :param tempdir: temporary directory
        :return:
        """
        infile = ds.GetDescription()

        import tempfile

        # record source projection into file, so these files are used as arguments for gdal warp later
        (fileh, src_prf) = tempfile.mkstemp('.prf', dir=tempdir)
        fileobj = os.fdopen(fileh, 'w')
        srcproj = ds.GetProjection()
        fileobj.write(srcproj)
        fileobj.close()

        # record target projection into file, so these files are used as arguments for gdal warp later
        (fileh, dest_prf) = tempfile.mkstemp('.prf', dir=tempdir)
        fileobj = os.fdopen(fileh, 'w')
        fileobj.write(self.reference_rast_grid.projection)
        fileobj.close()


        # use vrt as temporary file?
        ext = '.vrt'
        driver_name = 'VRT'

        # Now assemble command for gdalwarp
        # see http://www.gdal.org/gdalwarp.html for more info

        (fileh, temp_image) = tempfile.mkstemp(ext, dir=tempdir)
        os.close(fileh)

        cmd_list = [constants.GDALWARP, '-s_srs', src_prf, '-t_srs', dest_prf,
                    '-te',repr(working_region.min_x), # set georeferenced extents of output file to be created
                          repr(working_region.min_y), # (in target SRS by default, or in the SRS specified with -te_srs)
                          repr(working_region.max_x),
                          repr(working_region.max_y),
                    '-tr',repr(working_region.x_res), # set output file resolution (in target georeferenced units)
                          repr(working_region.y_res),
                    '-of',driver_name, # driver name to create file (VRT in this case)
                    '-r',resamplemethod] # reample method 'near' by default


        # optionnaly sets -srcnodata', '-dstnodata' if null_val_list provides them
        null_options = self.get_warp_null_options(null_val_list)
        if null_options is not None:
            cmd_list.extend(null_options)

        # finally, add input file
        cmd_list.append(infile)

        # and output file
        cmd_list.append(temp_image)

        # remember temporary files, so they are removed later
        self.temporary_files.append(temp_image)
        self.temporary_files.append(src_prf)
        self.temporary_files.append(dest_prf)

        # check that everything went smoothly during warping by return code
        ret_code = subprocess.call(cmd_list,
                                     stdout=self.logstream,
                                     stderr=self.logstream)

        # alarm, gdalwarp failed!
        if ret_code != 0:
            msg = 'gdalwarp exited with error code: {}'.format(ret_code)
            raise errors.GdalwarpError(msg)

        # check that resampled dataset opens just fine
        newds = gdal.Open(temp_image)
        if newds is None:
            msg = 'resampling failed, unable to open temporary image {}'.format(temp_image)
            raise errors.OpenImageError(msg)

        return newds

    def resample_all_to_reference(self, footprint, resamplemethodlist, tempdir='.'):
        """
        Resamples all files in container to same grid
        :param footprint: FP_UNION, FP_INTERSECTION or FP_REFERENCE
        :param resamplemethodlist: list of resample methods for each file
        :param tempdir: directory to store .vrt files into
        :return:
        """

        # find working region based on footprint
        # working region - raster grid defn
        working_region = self.find_working_region(footprint)

        # for each dataset in container
        for count in range(len(self.dataset_list)):

            # get raster grid
            rast_grid = self.rastgrid_list[count]

            # check if resampling is needed
            no_resample_needed = True
            if not self.reference_rast_grid.equal_resolution(rast_grid):
                self.logstream.write("Images Pixel size not equal {.20f} {.20f} {.20f} {.20f}\n".format(
                    self.reference_rast_grid.x_res, rast_grid.x_res, self.reference_rast_grid.y_res, rast_grid.y_res))
                no_resample_needed = False
            elif not self.reference_rast_grid.equal_projection(rast_grid):
                self.logstream.write("Images Coordinate system not equal {} {}\n".format(
                    self.reference_rast_grid.projection, rast_grid.projection))
                no_resample_needed = False
            elif not self.reference_rast_grid.aligned_with(rast_grid):
                self.logstream.write("Images grids aren't equal\n")
                no_resample_needed = False

            # yes, resampling is needed
            if not no_resample_needed:
                ds = self.dataset_list[count]
                null_vals = self.null_val_list[count]

                # resample by specific method
                resamplemethod = resamplemethodlist[count]

                # resampled dataset
                newds = self.resample_to_reference(ds, null_vals, working_region, resamplemethod, tempdir)

                # replace old adataset with new one
                self.dataset_list[count] = newds

                # replace old raster grid with new one
                self.rastgrid_list[count] = working_region

                del ds

    def is_all_match(self):
        """
        Check weather all dataset geocoordinates match. else, resampling is needed
        :return: true if all datasets have same working region
        """
        match = True
        for rast_grid in self.rastgrid_list:
            if not self.reference_rast_grid.equal_resolution(rast_grid):
                match = False
                break
            elif not self.reference_rast_grid.equal_projection(rast_grid):
                match = False
                break
            elif not self.reference_rast_grid.aligned_with(rast_grid):
                match = False
                break

        return match

    @staticmethod
    def get_warp_null_options(null_val_list):
        """
        internal function.
        Create options for gdalwarp based on null values of input and output datasets
        :param null_val_list:
        :return: list with -srcnodata and -dstnodata parameters, or None
        """
        have_none = False
        for null_val in null_val_list:
            if null_val is None:
                have_none = True
        if have_none:
            option_list = None
        else:
            null_val_str_list = [str(n) for n in null_val_list]
            all_nulls = ' '.join(null_val_str_list)
            option_list = ['-srcnodata', all_nulls, '-dstnodata', all_nulls]
        return option_list

    @staticmethod
    def make_rast_grid_from_dataset(ds):
        """
        Creates raster grid from dataset
        :param ds:
        :return:
        """
        proj = ds.GetProjection()
        geotrans = ds.GetGeoTransform()
        (nrows, ncols) = (ds.RasterYSize, ds.RasterXSize)
        rast_grid = rastergrid.RasterGrid(projection=proj, geotransform=geotrans, row_number=nrows,
                                          column_number=ncols)
        return rast_grid

    def find_working_region(self, footprint):
        """
        find common raster grid based on footprint option
        :param footprint: enum FP_UNION, FP_INTERSECTION, FP_REFERENCE
        :return: rester grid
        """

        combined_grid = rastergrid.find_common_region(self.rastgrid_list,
                                                      self.reference_rast_grid, footprint=footprint)
        return combined_grid
