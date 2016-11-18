# Load an ImageStack TMP file, an up to 4-dimensional array of integer or
# floating point data.
# Format specification is here:
# https://github.com/abadams/ImageStack/blob/master/src/FileTMP.cpp
#
# The output is permutted into the numpy C-contiguous convention: (frames, height, width, channels).
#   channels changes the fastest in memory, frames the slowest.
import numpy as np

def loadTMP(filename):
    with open(filename, 'rb') as f:
        buffer = f.read()

    type_code_list = \
    [
        np.dtype(np.float32),
        np.dtype(np.float64),
        np.dtype(np.uint8),
        np.dtype(np.int8),
        np.dtype(np.uint16),
        np.dtype(np.int16),
        np.dtype(np.uint32),
        np.dtype(np.int32),
        np.dtype(np.uint64),
        np.dtype(np.int64),
    ]

    offset = 0
    tmp_dims = np.frombuffer(buffer, count=4, dtype=np.int32)
    offset += tmp_dims.nbytes
    type_code = np.frombuffer(buffer, count=1, offset=offset, dtype=np.int32)
    offset += type_code.nbytes
    data_type = type_code_list[int(type_code)]

    # The TMP file writes the dimensions as:
    #   [width height frames channels]
    # Rewrite it in the numpy order:
    #   (channels, frames, height, width)
    a2_dims = np.flipud(tmp_dims)

    # Buffer stores the data as: (channels, frames, height, width).
    a2 = np.frombuffer(buffer, offset=offset, dtype=data_type)
    a2.shape = a2_dims

    # Permute the dimensions so that the data is returned in the numpy C-contiguous convention:
    #   (frames, height, width, channels).
    return a2.transpose((1, 2, 3, 0))

