import numpy as np
from loadTMP import loadTMP
from saveTMP import saveTMP

# Shift a uniformly distributed array to an integer range.
def _shift_to_int(float_array, lo, hi, dtype):
    return (lo + (hi - lo) * float_array).astype(dtype)

# Test TMP file reading and writing by saving to the given filename,
# reading it back, and checking for equality.
def testTMP(temp_filename):

    # Some decent size that's not too big.
    # Pick odd sizes on purpose.
    height = 767
    width = 1021
    channels = 3
    frames = 7
    int_lo = 37
    int_hi = 124    

    shape = (height, width, channels, frames)

    print('Testing float32...', end='')
    a = np.random.random_sample(shape).astype(np.float32)
    saveTMP(a, temp_filename)
    b = loadTMP(temp_filename)

    if np.array_equal(a, b):
        print('passed.')
    else:
        print('FAILED.')
        return False

    print('Testing double...', end='')
    a = np.random.random_sample(shape).astype(np.float64)
    saveTMP(a, temp_filename)
    b = loadTMP(temp_filename)

    if np.array_equal(a, b):
        print('passed.')
    else:
        print('FAILED.')
        return False

    print('Testing uint8...', end='')
    a = _shift_to_int(np.random.random_sample(shape),
                      int_lo, int_hi, np.uint8)
    saveTMP(a, temp_filename)
    b = loadTMP(temp_filename)

    if np.array_equal(a, b):
        print('passed.')
    else:
        print('FAILED.')
        return False

    print('Testing int8...', end='')
    a = _shift_to_int(np.random.random_sample(shape),
                      int_lo, int_hi, np.int8)
    saveTMP(a, temp_filename)
    b = loadTMP(temp_filename)

    if np.array_equal(a, b):
        print('passed.')
    else:
        print('FAILED.')
        return False

    print('Testing uint16...', end='')
    a = _shift_to_int(np.random.random_sample(shape),
                      int_lo, int_hi, np.uint16)
    saveTMP(a, temp_filename)
    b = loadTMP(temp_filename)

    if np.array_equal(a, b):
        print('passed.')
    else:
        print('FAILED.')
        return False

    print('Testing int16...', end='')
    a = _shift_to_int(np.random.random_sample(shape),
                      int_lo, int_hi, np.int16)
    saveTMP(a, temp_filename)
    b = loadTMP(temp_filename)

    if np.array_equal(a, b):
        print('passed.')
    else:
        print('FAILED.')
        return False

    print('Testing uint32...', end='')
    a = _shift_to_int(np.random.random_sample(shape),
                      int_lo, int_hi, np.uint32)
    saveTMP(a, temp_filename)
    b = loadTMP(temp_filename)

    if np.array_equal(a, b):
        print('passed.')
    else:
        print('FAILED.')
        return False

    print('Testing int32...', end='')
    a = _shift_to_int(np.random.random_sample(shape),
                      int_lo, int_hi, np.int32)
    saveTMP(a, temp_filename)
    b = loadTMP(temp_filename)

    if np.array_equal(a, b):
        print('passed.')
    else:
        print('FAILED.')
        return False

    print('Testing uint64...', end='')
    a = _shift_to_int(np.random.random_sample(shape),
                      int_lo, int_hi, np.uint64)
    saveTMP(a, temp_filename)
    b = loadTMP(temp_filename)

    if np.array_equal(a, b):
        print('passed.')
    else:
        print('FAILED.')
        return False

    print('Testing int64...', end='')
    a = _shift_to_int(np.random.random_sample(shape),
                      int_lo, int_hi, np.int64)
    saveTMP(a, temp_filename)
    b = loadTMP(temp_filename)

    if np.array_equal(a, b):
        print('passed.')
    else:
        print('FAILED.')
        return False

    print('All tests PASSED.')
    return True
