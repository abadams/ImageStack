# Write an up to 4-dimensional numpy.array 'a' to 'filename' in the ImageStack TMP format.
# Format specification is here:
# https://github.com/abadams/ImageStack/blob/master/src/FileTMP.cpp
#
# We assume a is in numpy's C-continguous convention: (frames, height, width, channels).
#   channels changes the fastest in memory, frames the slowest.
#   It does not have to be stored that way, we just assume those are the dimensions.
#
# a will be written to disk in the TMP file convention: (channels frames height width).
#   width changes the fastest in memory, channels the slowest.
#   Note that the TMP file writes dimensions in the opposite order from numpy:
#   [width, height, frames, channels].
#
import numpy as np

def saveTMP(a, filename):
    if a.ndim > 4:
        raise ValueError('TMP only supports up to 4 dimensions.')

    type_code_map = \
    {
        np.dtype(np.float32): np.int32(0),
        np.dtype(np.float64): np.int32(1),
        np.dtype(np.uint8): np.int32(2),
        np.dtype(np.int8): np.int32(3),
        np.dtype(np.uint16): np.int32(4),
        np.dtype(np.int16): np.int32(5),
        np.dtype(np.uint32): np.int32(6),
        np.dtype(np.int32): np.int32(7),
        np.dtype(np.uint64): np.int32(8),
        np.dtype(np.int64): np.int32(9),
    }

    type_code = type_code_map.get(a.dtype, np.int32(-1))
    if type_code == -1:
        raise ValueError(
            'TMP only supports float{32,64} and {u}int{8,16,32,64} numpy.arrays.')

    # Permute the axes to [width height frames channels].
    # First pad the dimensions
    a2 = a.copy()
    # f.atleast_3d is nice, but numpy doesn't come with a 4D version.
    # I also want to add dimensions to the front, not the end.
    a2.shape = (1,) * (4 - a2.ndim) + a2.shape
    a2 = a2.transpose((3, 0, 1, 2))

    # TMP wants dimensions written in the opposite order.
    tmp_dims = np.flipud(np.array(a2.shape, dtype=np.int32))

    with open(filename, 'wb') as f:
        f.write(tmp_dims.tobytes())
        f.write(type_code.tobytes())
        f.write(a2.tobytes())
