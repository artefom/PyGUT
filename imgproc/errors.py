# Artyom Fomenko BSE-142


class ImgprocError(Exception): pass


class ImageSizeError(ImgprocError):
    "demanded block is out of range"


class GdalwarpError(ImgprocError):
    "gdalwarp not found"


class OpenFileError(ImgprocError):
    "Could not open resulting or input file"


class OpenImageError(OpenFileError):
    "GDAL failed to read image"


class InvalidArgumentsError(ImgprocError):
    "Function got incorrect arguments"


class ResampleRequiedError(ImgprocError):
    "Image's coordinates don't match. resample is requied"


class InvalidKey(ImgprocError):
    "Unexpected keys encountered"


class InvalidType(ImgprocError):
    "Failed to convert type"


class ListLengthError(ImgprocError):
    "Lists length doesn't match"


class InvalidArrayShape(ImgprocError):
    "Invalid array shape encountered"


class InvalidIntersection(ImgprocError):
    "Images do not intersect"
