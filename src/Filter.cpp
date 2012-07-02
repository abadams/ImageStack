#include "main.h"
#include "Filter.h"
#include "Convolve.h"
#include "Color.h"
#include "Geometry.h"
#include "Arithmetic.h"
#include "Statistics.h"
#include "header.h"

void GaussianBlur::help() {
    pprintf("-gaussianblur takes a floating point width, height, and frames, and"
            " performs a gaussian blur with those standard deviations. The blur is"
            " performed out to three standard deviations. If given only two"
            " arguments, it performs a blur in x and y only. If given one argument,"
            " it performs the blur in x and y with filter width the same as"
            " height.\n"
            "\n"
            "Usage: ImageStack -load in.jpg -gaussianblur 5 -save blurry.jpg\n\n");
}

bool GaussianBlur::test() {
    Image impulse(21, 21, 21, 3);
    impulse(10, 10, 10, 0) = 1;
    impulse(10, 10, 10, 1) = 2;
    impulse(10, 10, 10, 2) = 3;
    Image blurry = GaussianBlur::apply(impulse, 0.7, 0.8, 0.5);
    float ratio = blurry(10, 10, 10, 0);
    for (int t = 0; t < 21; t++) {
        float ft = (t - 10.0f)/0.5f;
        for (int y = 0; y < 21; y++) {
            float fy = (y - 10.0f)/0.8f;
            for (int x = 0; x < 21; x++) {
                float fx = (x - 10.0f)/0.7f;
                float correct = expf(-0.5f*(fx*fx + fy*fy + ft*ft))*ratio;
                if (!nearlyEqual(blurry(x, y, t, 0), correct*1)) return false;
                if (!nearlyEqual(blurry(x, y, t, 1), correct*2)) return false;
                if (!nearlyEqual(blurry(x, y, t, 2), correct*3)) return false;
            }
        }
    }
    return nearlyEqual(Stats(blurry).sum(), 6);
}

void GaussianBlur::parse(vector<string> args) {
    float frames = 0, width = 0, height = 0;
    if (args.size() == 1) {
        width = height = readFloat(args[0]);
    } else if (args.size() == 2) {
        width = readFloat(args[0]);
        height = readFloat(args[1]);
    } else if (args.size() == 3) {
        width  = readFloat(args[0]);
        height = readFloat(args[1]);
        frames = readFloat(args[2]);
    } else {
        panic("-gaussianblur takes one, two, or three arguments\n");
    }

    Image im = apply(stack(0), width, height, frames);
    pop();
    push(im);
}

Image GaussianBlur::apply(Image im, float filterWidth, float filterHeight, float filterFrames) {
    Image out(im);

    if (filterWidth != 0) {
        // make the width filter
        int size = (int)(filterWidth * 6 + 1) | 1;
        // even tiny filters should do something, otherwise we
        // wouldn't have called this function.
        if (size == 1) { size = 3; }
        int radius = size / 2;
        Image filter(size, 1, 1, 1);
        float sum = 0;
        for (int i = 0; i < size; i++) {
            float diff = (i-radius)/filterWidth;
            float value = expf(-diff * diff / 2);
            filter(i, 0, 0, 0) = value;
            sum += value;
        }

        for (int i = 0; i < size; i++) {
            filter(i, 0, 0, 0) /= sum;
        }

        out = Convolve::apply(out, filter, Convolve::Homogeneous);
    }

    if (filterHeight != 0) {
        // make the height filter
        int size = (int)(filterHeight * 6 + 1) | 1;
        // even tiny filters should do something, otherwise we
        // wouldn't have called this function.
        if (size == 1) { size = 3; }
        int radius = size / 2;
        Image filter(1, size, 1, 1);
        float sum = 0;
        for (int i = 0; i < size; i++) {
            float diff = (i-radius)/filterHeight;
            float value = expf(-diff * diff / 2);
            filter(0, i, 0, 0) = value;
            sum += value;
        }

        for (int i = 0; i < size; i++) {
            filter(0, i, 0, 0) /= sum;
        }

        out = Convolve::apply(out, filter, Convolve::Homogeneous);
    }

    if (filterFrames != 0) {
        // make the frames filter
        int size = (int)(filterFrames * 6 + 1) | 1;
        // even tiny filters should do something, otherwise we
        // wouldn't have called this function.
        if (size == 1) { size = 3; }
        int radius = size / 2;
        Image filter(1, 1, size, 1);
        float sum = 0;
        for (int i = 0; i < size; i++) {
            float diff = (i-radius)/filterFrames;
            float value = expf(-diff * diff / 2);
            filter(0, 0, i, 0) = value;
            sum += value;
        }

        for (int i = 0; i < size; i++) {
            filter(0, 0, i, 0) /= sum;
        }

        out = Convolve::apply(out, filter, Convolve::Homogeneous);
    }

    return out;
}

// This blur implementation was contributed by Tyler Mullen as a
// CS448F project. A competition was held, and this method was found
// to be much faster than other IIRs, filtering by resampling,
// iterated rect filters, and polynomial integral images. The method
// was modified by Andrew Adams to be more ImageStacky (i.e. use
// structures more idiomatic to ImageStack), to work for larger sized
// blurs, and to cover more unusual cases.

void FastBlur::help() {
    pprintf("-fastblur takes a floating point width, height, and frames, and"
            " performs a fast approximate gaussian blur with those standard"
            " deviations using the IIR method of van Vliet et al. If given only two"
            " arguments, it performs a blur in x and y only. If given one argument,"
            " it performs the blur in x and y with filter width the same as"
            " height.\n"
            "\n"
            "Usage: ImageStack -load in.jpg -fastblur 5 -save blurry.jpg\n\n");
}

bool FastBlur::test() {
    Image a(100, 80, 40, 4);
    Noise::apply(a, 0, 1);
    Image b = GaussianBlur::apply(a, 2.3, 1.3, 1.2);
    FastBlur::apply(a, 2.3, 1.3, 1.2);
    return nearlyEqual(a, b);
}

void FastBlur::parse(vector<string> args) {
    float frames = 0, width = 0, height = 0;
    if (args.size() == 1) {
        width = height = readFloat(args[0]);
    } else if (args.size() == 2) {
        width = readFloat(args[0]);
        height = readFloat(args[1]);
    } else if (args.size() == 3) {
        width  = readFloat(args[0]);
        height = readFloat(args[1]);
        frames = readFloat(args[2]);
    } else {
        panic("-fastblur takes one, two, or three arguments\n");
    }

    apply(stack(0), width, height, frames);
}

void FastBlur::apply(Image im, float filterWidth, float filterHeight, float filterFrames) {
    assert(filterFrames >= 0 &&
           filterWidth >= 0 &&
           filterHeight >= 0,
           "Filter sizes must be non-negative\n");

    // Prevent filtering in useless directions
    if (im.width == 1) { filterWidth = 0; }
    if (im.height == 1) { filterHeight = 0; }
    if (im.frames == 1) { filterFrames = 0; }

    //printf("%d %d %d\n", im.width, im.height, im.frames);

    // Filter in very narrow directions using the regular Gaussian, as
    // the IIR requires a few pixels to get going. If the Gaussian
    // blur is very narrow, also revert to the naive method, as IIR
    // won't work.
    if (filterWidth > 0 && (im.width < 16 || filterWidth < 0.5)) {
        Image blurry = GaussianBlur::apply(im, filterWidth, 0, 0);
        FastBlur::apply(blurry, 0, filterHeight, filterFrames);
        im.set(blurry);
        return;
    }

    if (filterHeight > 0 && (im.height < 16 || filterHeight < 0.5)) {
        Image blurry = GaussianBlur::apply(im, 0, filterHeight, 0);
        FastBlur::apply(blurry, filterWidth, 0, filterFrames);
        im.set(blurry);
        return;
    }

    if (filterFrames > 0 && (im.frames < 16 || filterFrames < 0.5)) {
        Image blurry = GaussianBlur::apply(im, 0, 0, filterFrames);
        FastBlur::apply(blurry, filterWidth, filterHeight, 0);
        im.set(blurry);
        return;
    }

    // Deal with very large filters by splitting into multiple smaller filters
    int xIterations = 1, yIterations = 1, tIterations = 1;
    while (filterWidth > 64) {
        filterWidth /= sqrtf(2);
        xIterations *= 2;
    }
    while (filterHeight > 64) {
        filterHeight /= sqrtf(2);
        yIterations *= 2;
    }
    while (filterFrames > 64) {
        filterFrames /= sqrtf(2);
        tIterations *= 2;
    }

    const int w = 16;

    // blur in x
    if (filterWidth > 0) {
        const int size = im.width + (int)(filterWidth*6);

        float c0, c1, c2, c3;
        calculateCoefficients(filterWidth, &c0, &c1, &c2, &c3);

        vector<float> scale(size);
        computeAttenuation(&scale[0], size, im.width, c0, c1, c2, c3, xIterations);

        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
#ifdef _OPENMP
                #pragma omp parallel for
#endif
                for (int y = 0; y < im.height; y += w) {
                    vector<float> chunk(size*w, 0);

                    // prepare 16 scanlines
                    for (int x = 0; x < im.width; x++) {
                        for (int i = 0; i < w && y+i < im.height; i++) {
                            chunk[x*w + i] = im(x, y+i, t, c);
                        }
                    }

                    // blur them
                    for (int i = 0; i < xIterations; i++) {
                        blurChunk(&chunk[0], size, c0, c1, c2, c3);
                        blurChunk(&chunk[0], size, c0, c1, c2, c3);
                    }

                    // read them back
                    for (int x = 0; x < im.width; x++) {
                        for (int i = 0; i < w && y+i < im.height; i++) {
                            im(x, y+i, t, c) = chunk[x*w + i] * scale[x];
                        }
                    }
                }
            }
        }
    }

    // blur in y
    if (filterHeight > 0) {
        const int size = im.height + (int)(filterHeight*6);

        float c0, c1, c2, c3;
        calculateCoefficients(filterHeight, &c0, &c1, &c2, &c3);

        vector<float> scale(size);
        computeAttenuation(&scale[0], size, im.height, c0, c1, c2, c3, yIterations);

        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
#ifdef _OPENMP
                #pragma omp parallel for
#endif
                for (int x = 0; x < im.width; x += w) {
                    vector<float> chunk(size*w, 0);

                    // prepare 16 columns
                    for (int y = 0; y < im.height; y++) {
                        for (int i = 0; i < w && x+i < im.width; i++) {
                            chunk[y*w + i] = im(x+i, y, t, c);
                        }
                    }

                    // blur them
                    for (int i = 0; i < yIterations; i++) {
                        blurChunk(&chunk[0], size, c0, c1, c2, c3);
                        blurChunk(&chunk[0], size, c0, c1, c2, c3);
                    }

                    // read them back
                    for (int y = 0; y < im.height; y++) {
                        for (int i = 0; i < w && x+i < im.width; i++) {
                            im(x+i, y, t, c) = chunk[y*w + i] * scale[y];
                        }
                    }
                }
            }
        }
    }

    // blur in t
    if (filterFrames > 0) {
        const int size = im.frames + (int)(filterFrames*6);

        float c0, c1, c2, c3;
        calculateCoefficients(filterFrames, &c0, &c1, &c2, &c3);

        vector<float> scale(size);
        computeAttenuation(&scale[0], size, im.frames, c0, c1, c2, c3, tIterations);

        for (int c = 0; c < im.channels; c++) {
            for (int y = 0; y < im.height; y++) {
#ifdef _OPENMP
                #pragma omp parallel for
#endif
                for (int x = 0; x < im.width; x += w) {
                    vector<float> chunk(size*w, 0);

                    // prepare 16 scanlines
                    for (int t = 0; t < im.frames; t++) {
                        for (int i = 0; i < w && x+i < im.width; i++) {
                            chunk[t*w + i] = im(x+i, y, t, c);
                        }
                    }

                    // blur them
                    for (int i = 0; i < tIterations; i++) {
                        blurChunk(&chunk[0], size, c0, c1, c2, c3);
                        blurChunk(&chunk[0], size, c0, c1, c2, c3);
                    }

                    // read them back
                    for (int t = 0; t < im.frames; t++) {
                        for (int i = 0; i < w && x+i < im.width; i++) {
                            im(x+i, y, t, c) = chunk[t*w + i] * scale[t];
                        }
                    }
                }
            }
        }
    }
}

void FastBlur::blurChunk(float *data, int size,
                         const float c0, const float c1,
                         const float c2, const float c3) {
    // filter a 16-wide chunk of data in place

    const int w = 16;

    // Warm up
    for (int x = 0; x < w; x++) {
        data[x] = c0*data[x];
        data[w + x] = (c0*data[w + x] +
                       c1*data[x]);
        data[2*w + x] = (c0*data[2*w + x] +
                         c1*data[w + x] +
                         c2*data[x]);
    }

    // Filter
    for (int x = 3*w; x < size*w; x++) {
        data[x] = (c0 * data[x] +
                   c1 * data[x-w] +
                   c2 * data[x-w*2] +
                   c3 * data[x-w*3]);
    }

    // Flip the data
    for (int y = 0; y < size/2; y++) {
        for (int x = 0; x < w; x++) {
            swap(data[y*w+x], data[(size-1-y)*w+x]);
        }
    }
}

void FastBlur::computeAttenuation(float *scale, int size, int width,
                                  const float c0, const float c1,
                                  const float c2, const float c3,
                                  int iterations) {
    // Initial value
    for (int x = 0; x < width; x++) {
        scale[x] = 1.0f;
    }
    for (int x = width; x < size; x++) {
        scale[x] = 0.0f;
    }

    for (int i = 0; i < iterations; i++) {
        // Forward pass
        scale[0] = c0 * scale[0];
        scale[1] = c0 * scale[1] + c1 * scale[0];
        scale[2] = c0 * scale[2] + c1 * scale[1] + c2 * scale[0];
        for (int x = 3; x < size; x++) {
            scale[x] = (c0 * scale[x] +
                        c1 * scale[x-1] +
                        c2 * scale[x-2] +
                        c3 * scale[x-3]);
        }

        // Backward pass
        scale[size-1] = c0 * scale[size-1];
        scale[size-2] = (c0 * scale[size-2] +
                         c1 * scale[size-1]);
        scale[size-3] = (c0 * scale[size-3] +
                         c1 * scale[size-2] +
                         c2 * scale[size-1]);
        for (int x = size-4; x >= 0; x--) {
            scale[x] = (c0 * scale[x] +
                        c1 * scale[x+1] +
                        c2 * scale[x+2] +
                        c3 * scale[x+3]);
        }
    }

    // Invert
    for (int x = 0; x < size; x++) {
        scale[x] = 1.0f/scale[x];
    }
}

void FastBlur::calculateCoefficients(float sigma, float *c0, float *c1, float *c2, float *c3) {
    // performs the necessary conversion between the sigma of a Gaussian blur
    // and the coefficients used in the IIR filter

    float q;

    assert(sigma >= 0.5, "To use IIR filtering, standard deviation of blur must be >= 0.5\n");

    if (sigma < 2.5) {
        q = (3.97156 - 4.14554*sqrtf(1 - 0.26891*sigma));
    } else {
        q = 0.98711*sigma - 0.96330;
    }

    float denom = 1.57825 + 2.44413*q + 1.4281*q*q + 0.422205*q*q*q;
    *c1 = (2.44413*q + 2.85619*q*q + 1.26661*q*q*q)/denom;
    *c2 = -(1.4281*q*q + 1.26661*q*q*q)/denom;
    *c3 = (0.422205*q*q*q)/denom;
    *c0 = 1 - (*c1 + *c2 + *c3);
}

void RectFilter::help() {
    pprintf("-rectfilter performs a iterated rectangular filter on the image. The"
            " four arguments are the filter width, height, frames, and the number of"
            " iterations. If three arguments are given, they are interpreted as"
            " frames, width, and height, and the number of iterations is assumed to"
            " be one. If two arguments are given they are taken as width and height,"
            " and frames is assumed to be one. If one argument is given it is taken"
            " as both width and height, with frames and iterations again assumed to"
            " be one.\n"
            "\n"
            "Usage: ImageStack -load in.jpg -rectfilter 1 10 10 -save out.jpg\n\n");
}

bool RectFilter::test() {
    Image impulse(21, 21, 21, 3);
    impulse(10, 10, 10, 0) = 1;
    impulse(10, 10, 10, 1) = 2;
    impulse(10, 10, 10, 2) = 3;
    RectFilter::apply(impulse, 3, 5, 7, 1);
    for (int t = 0; t < 21; t++) {
        for (int y = 0; y < 21; y++) {
            for (int x = 0; x < 21; x++) {
                bool inside = (x > 8 && x < 12 &&
                               y > 7 && y < 13 &&
                               t > 6 && t < 14);
                float correct = inside ? 1.0f/(3*5*7) : 0.0f;
                if (!nearlyEqual(impulse(x, y, t, 0), correct*1)) return false;
                if (!nearlyEqual(impulse(x, y, t, 1), correct*2)) return false;
                if (!nearlyEqual(impulse(x, y, t, 2), correct*3)) return false;
            }
        }
    }
    return true;
}

void RectFilter::parse(vector<string> args) {
    int iterations = 1, frames = 1, width = 1, height = 1;
    if (args.size() == 1) {
        width = height = readInt(args[0]);
    } else if (args.size() == 2) {
        width = readInt(args[0]);
        height = readInt(args[1]);
    } else if (args.size() == 3) {
        width = readInt(args[0]);
        height = readInt(args[1]);
        frames = readInt(args[2]);
    } else if (args.size() == 4) {
        width = readInt(args[0]);
        height = readInt(args[1]);
        frames = readInt(args[2]);
        iterations = readInt(args[3]);
    } else {
        panic("-rectfilter takes four or fewer arguments\n");
    }

    apply(stack(0), width, height, frames, iterations);
}

void RectFilter::apply(Image im, int filterWidth, int filterHeight, int filterFrames, int iterations) {
    assert(filterFrames & filterWidth & filterHeight & 1, "filter shape must be odd\n");
    assert(iterations >= 1, "iterations must be at least one\n");

    if (filterFrames != 1) blurT(im, filterFrames, iterations);
    if (filterWidth  != 1) blurX(im, filterWidth, iterations);
    if (filterHeight != 1) blurY(im, filterHeight, iterations);
}

void RectFilter::blurXCompletely(Image im) {
    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                // compute the average for this scanline
                double average = 0;
                for (int x = 0; x < im.width; x++) {
                    average += im(x, y, t, c);
                }
                average /= im.width;
                for (int x = 0; x < im.width; x++) {
                    im(x, y, t, c) = average;
                }
            }
        }
    }
}


void RectFilter::blurX(Image im, int width, int iterations) {
    if (width <= 1) { return; }
    if (im.width == 1) { return; }

    // special case where the radius is large enough that the image is totally uniformly blurred
    if (im.width <= width/2) {
        blurXCompletely(im);
        return;
    }

    int radius = width/2;
    vector<float> buffer(width);

    vector<float> multiplier(width);
    for (int i = 0; i < width; i++) {
        multiplier[i] = 1.0f/width;
    }

    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                for (int i = 0; i < iterations; i++) {
                    // keep a circular buffer of everything currently inside the kernel
                    // also maintain the sum of this buffer

                    double sum = 0;
                    int bufferIndex = 0;
                    int bufferEntries = 0;

                    // initialize the buffer
                    for (int j = 0; j <= radius; j++) {
                        buffer[j] = 0;
                    }
                    for (int j = radius+1; j < width; j++) {
                        buffer[j] = im(j-radius, y, t, c);
                        sum += buffer[j];
                        bufferEntries++;
                    }

                    double mult = 1.0/bufferEntries;

                    // non boundary cases
                    for (int x = 0; x < im.width-radius-1; x++) {
                        // assign the average to the current position
                        im(x, y, t, c) = (float)(sum * mult);

                        // swap out the buffer element, updating the sum
                        float newVal = im(x+radius+1, y, t, c);
                        sum += newVal - buffer[bufferIndex];
                        buffer[bufferIndex] = newVal;
                        bufferIndex++;
                        if (bufferIndex == width) { bufferIndex = 0; }

                        if (bufferEntries < width) {
                            bufferEntries++;
                            mult = 1.0/bufferEntries;
                        }
                    }

                    // boundary cases
                    for (int x = im.width-radius-1; x < im.width; x++) {
                        // assign the average to the current position
                        im(x, y, t, c) = (float)(sum * mult);

                        // swap out the buffer element, updating the sum
                        sum -= buffer[bufferIndex];
                        //buffer[bufferIndex] = 0;
                        bufferIndex++;
                        if (bufferIndex == width) { bufferIndex = 0; }

                        bufferEntries--;
                        mult = 1.0/bufferEntries;
                    }
                }
            }
        }
    }

}

void RectFilter::blurY(Image im, int width, int iterations) {
    if (width <= 1) { return; }
    if (im.height == 1) { return; }

    // pull out strips of columns and blur them
    Image chunk(im.height, 8, 1, 1);

    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames; t++) {
            for (int x = 0; x < im.width; x += chunk.height) {
                int size = chunk.height;
                if (x + chunk.height >= im.width) { size = im.width-x; }

                // read into the chunk in a transposed fashion
                for (int y = 0; y < im.height; y++) {
                    for (int j = 0; j < size; j++) {
                        chunk(y, j) = im(x+j, y, t, c);
                    }
                }

                // blur the chunk
                blurX(chunk, width, iterations);

                // read back from the chunk
                for (int y = 0; y < im.height; y++) {
                    for (int j = 0; j < size; j++) {
                        im(x+j, y, t, c) = chunk(y, j);
                    }
                }
            }
        }
    }
}

void RectFilter::blurT(Image im, int width, int iterations) {
    if (width <= 1) { return; }
    if (im.frames == 1) { return; }

    // pull out strips across frames from rows and blur them
    Image chunk(im.frames, 8, 1, 1);

    for (int c = 0; c < im.channels; c++) {
        for (int y = 0; y < im.height; y++) {
            for (int x = 0; x < im.width; x += chunk.height) {
                int size = chunk.height;
                if (x + chunk.height >= im.width) { size = im.width-x; }

                // read into the chunk in a transposed fashion
                for (int t = 0; t < im.frames; t++) {
                    for (int j = 0; j < size; j++) {
                        chunk(t, j) = im(x+j, y, t, c);
                    }
                }

                // blur the chunk
                blurX(chunk, width, iterations);

                // read back from the chunk
                for (int t = 0; t < im.frames; t++) {
                    for (int j = 0; j < size; j++) {
                        im(x+j, y, t, c) = chunk(t, j);
                    }
                }
            }
        }
    }
}

void LanczosBlur::help() {
    pprintf("-lanczosblur convolves the current image by a three lobed lanczos"
            " filter. A lanczos filter is a kind of windowed sinc. The three"
            " arguments are filter width, height, and frames. If two arguments are"
            " given, frames is assumed to be one. If one argument is given, it is"
            " interpreted as both width and height.\n"
            "\n"
            "Usage: ImageStack -load big.jpg -lanczosblur 2 -subsample 2 2 0 0 -save small.jpg\n\n");

}

bool LanczosBlur::test() {
    Image impulse(21, 21, 21, 3);
    impulse(10, 10, 10, 0) = 1;
    impulse(10, 10, 10, 1) = 2;
    impulse(10, 10, 10, 2) = 3;
    Image blurry = LanczosBlur::apply(impulse, 1.7, 1.8, 1.5);
    float ratio = blurry(10, 10, 10, 0);
    for (int t = 0; t < 21; t++) {
        float ft = (t - 10.0f)/1.5f;
        ft = (ft == 0) ? 1 : (sinc(ft) * sinc(ft/3));
        for (int y = 0; y < 21; y++) {
            float fy = (y - 10.0f)/1.8f;
            fy = (fy == 0) ? 1 : (sinc(fy) * sinc(fy/3));
            for (int x = 0; x < 21; x++) {
                float fx = (x - 10.0f)/1.7f;
                fx = (fx == 0) ? 1 : (sinc(fx) * sinc(fx/3));
                float correct = fx*fy*ft*ratio;
                if (!nearlyEqual(blurry(x, y, t, 0), correct*1)) return false;
                if (!nearlyEqual(blurry(x, y, t, 1), correct*2)) return false;
                if (!nearlyEqual(blurry(x, y, t, 2), correct*3)) return false;
            }
        }
    }
    return nearlyEqual(Stats(blurry).sum(), 6);
}

void LanczosBlur::parse(vector<string> args) {
    float frames = 0, width = 0, height = 0;
    if (args.size() == 1) {
        width = height = readFloat(args[0]);
    } else if (args.size() == 2) {
        width = readFloat(args[0]);
        height = readFloat(args[1]);
    } else if (args.size() == 3) {
        width  = readFloat(args[0]);
        height = readFloat(args[1]);
        frames = readFloat(args[2]);
    } else {
        panic("-lanczosblur takes one, two, or three arguments\n");
    }

    Image im = apply(stack(0), width, height, frames);
    pop();
    push(im);
}

Image LanczosBlur::apply(Image im, float filterWidth, float filterHeight, float filterFrames) {
    Image out(im);

    if (filterFrames != 0) {
        // make the frames filter
        int size = (int)(filterFrames * 6 + 1) | 1;
        int radius = size / 2;
        Image filter(1, 1, size, 1);
        float sum = 0;
        for (int i = 0; i < size; i++) {
            float value = lanczos_3((i-radius) / filterFrames);
            filter(0, 0, i, 0) = value;
            sum += value;
        }

        for (int i = 0; i < size; i++) {
            filter(0, 0, i, 0) /= sum;
        }

        out = Convolve::apply(out, filter, Convolve::Homogeneous);
    }

    if (filterWidth != 0) {
        // make the width filter
        int size = (int)(filterWidth * 6 + 1) | 1;
        int radius = size / 2;
        Image filter(size, 1, 1, 1);
        float sum = 0;
        for (int i = 0; i < size; i++) {
            float value = lanczos_3((i-radius) / filterWidth);
            filter(i, 0, 0, 0) = value;
            sum += value;
        }

        for (int i = 0; i < size; i++) {
            filter(i, 0, 0, 0) /= sum;
        }

        out = Convolve::apply(out, filter, Convolve::Homogeneous);
    }

    if (filterHeight != 0) {
        // make the height filter
        int size = (int)(filterHeight * 6 + 1) | 1;
        int radius = size / 2;
        Image filter(1, size, 1, 1);
        float sum = 0;
        for (int i = 0; i < size; i++) {
            float value = lanczos_3((i-radius) / filterHeight);
            filter(0, i, 0, 0) = value;
            sum += value;
        }

        for (int i = 0; i < size; i++) {
            filter(0, i, 0, 0) /= sum;
        }

        out = Convolve::apply(out, filter, Convolve::Homogeneous);
    }

    return out;

}


void MinFilter::help() {
    pprintf("-minfilter applies a min filter with square support. The sole argument "
            "is the pixel radius of the filter. For circular support, see "
            "-percentilefilter.\n"
            "\n"
            "Usage: ImageStack -load input.jpg -minfilter 10 -save output.jpg\n");
}

bool MinFilter::test() {
    // check it's less than the input, but touches the input in some places
    Image a(123, 234, 2, 3);
    Noise::apply(a, 0, 1);
    Image b = a.copy();
    MinFilter::apply(b, 3);
    a -= b;
    return nearlyEqual(Stats(a).minimum(), 0);
}

void MinFilter::parse(vector<string> args) {
    assert(args.size() == 1, "-minfilter takes on argument\n");
    int radius = readInt(args[0]);
    assert(radius > -1, "radius must be positive");
    apply(stack(0), radius);
}

void MinFilter::apply(Image im, int radius) {
    // Make a heap with (2*radius + 1) leaves. Unlike a regular heap,
    // each internal node is a _copy_ of the smaller child. The leaf
    // nodes act as a circular buffer. Every time we introduce a new
    // pixel (and evict an old one), we update all of its parents up
    // to the root.

    vector<float> heap(4*radius+1);

    for (int t = 0; t < im.frames; t++) {
        for (int y = 0; y < im.height; y++) {
            for (int c = 0; c < im.channels; c++) {
                // Initialize the heap to contain all inf
                std::fill(heap.begin(), heap.end(), INF);
                size_t pos = 2*radius;
                for (int x = 0; x < im.width + radius; x++) {
                    // Get the next input
                    float val;
                    if (x < im.width) val = im(x, y, t, c);
                    else val = INF;
                    // Stuff it in the heap
                    heap[pos] = val;
                    // Update parents
                    size_t p = pos;
                    do {
                        p--; p>>=1;
                        heap[p] = min(heap[2*p+1], heap[2*p+2]);
                    } while (p);

                    // Maybe write out min
                    if (x-radius > 0)
                        im(x-radius, y, t, c) = heap[0];
                    // Update position in circular buffer
                    pos++;
                    if (pos == heap.size()) pos = 2*radius;
                }
            }
        }

        for (int x = 0; x < im.width; x++) {
            for (int c = 0; c < im.channels; c++) {
                // Initialize the heap to contain all inf
                std::fill(heap.begin(), heap.end(), INF);
                size_t pos = 2*radius;
                for (int y = 0; y < im.height + radius; y++) {
                    float val;
                    if (y < im.height) val = im(x, y, t, c);
                    else val = INF;
                    // stuff it in the heap
                    heap[pos] = val;
                    // update parents
                    size_t p = pos;
                    do {
                        p--; p>>=1;
                        heap[p] = min(heap[2*p+1], heap[2*p+2]);
                    } while (p);
                    // write out min
                    if (y-radius > 0)
                        im(x, y-radius, t, c) = heap[0];
                    // update position in circular buffer
                    pos++;
                    if (pos == heap.size()) pos = 2*radius;
                }
            }
        }
    }

}

void MaxFilter::help() {
    pprintf("-maxfilter applies a max filter with square support. The sole argument "
            "is the pixel radius of the filter. For circular support, see "
            "-percentilefilter.\n"
            "\n"
            "Usage: ImageStack -load input.jpg -maxfilter 10 -save output.jpg\n");
}

#include "File.h"

bool MaxFilter::test() {
    Image a(123, 234, 2, 3);
    Noise::apply(a, 0, 1);
    Image b = a.copy();
    MaxFilter::apply(b, 3);
    return nearlyEqual(Stats(b-a).minimum(), 0);
}

void MaxFilter::parse(vector<string> args) {
    assert(args.size() == 1, "-maxfilter takes on argument\n");
    int radius = readInt(args[0]);
    assert(radius > -1, "radius must be positive");
    apply(stack(0), radius);
}

void MaxFilter::apply(Image im, int radius) {
    // Make a heap with (2*radius + 1) leaves. Unlike a regular heap,
    // each internal node is a _copy_ of the smaller child. The leaf
    // nodes act as a circular buffer. Every time we introduce a new
    // pixel (and evict an old one), we update all of its parents up
    // to the root.

    vector<float> heap(4*radius+1);

    for (int t = 0; t < im.frames; t++) {
        for (int y = 0; y < im.height; y++) {
            for (int c = 0; c < im.channels; c++) {
                // Initialize the heap to contain all inf
                std::fill(heap.begin(), heap.end(), -INF);
                size_t pos = 2*radius;
                for (int x = 0; x < im.width + radius; x++) {
                    // Get the next input
                    float val;
                    if (x < im.width) val = im(x, y, t, c);
                    else val = -INF;
                    // Stuff it in the heap
                    heap[pos] = val;
                    // Update parents
                    size_t p = pos;
                    do {
                        p--; p>>=1;
                        heap[p] = max(heap[2*p+1], heap[2*p+2]);
                    } while (p);

                    // Maybe write out max
                    if (x-radius > 0)
                        im(x-radius, y, t, c) = heap[0];
                    // Update position in circular buffer
                    pos++;
                    if (pos == heap.size()) pos = 2*radius;
                }
            }
        }

        for (int x = 0; x < im.width; x++) {
            for (int c = 0; c < im.channels; c++) {
                // Initialize the heap to contain all inf
                std::fill(heap.begin(), heap.end(), -INF);
                size_t pos = 2*radius;
                for (int y = 0; y < im.height + radius; y++) {
                    float val;
                    if (y < im.height) val = im(x, y, t, c);
                    else val = -INF;
                    // Stuff it in the heap
                    heap[pos] = val;
                    // Update parents
                    size_t p = pos;
                    do {
                        p--; p>>=1;
                        heap[p] = max(heap[2*p+1], heap[2*p+2]);
                    } while (p);
                    // write out max
                    if (y-radius > 0)
                        im(x, y-radius, t, c) = heap[0];
                    // update position in circular buffer
                    pos++;
                    if (pos == heap.size()) pos = 2*radius;
                }
            }
        }
    }

}


void MedianFilter::help() {
    pprintf("-medianfilter applies a median filter with a circular support. The "
            "sole argument is the pixel radius of the filter.\n"
            "\n"
            "Usage: ImageStack -load input.jpg -medianfilter 10 -save output.jpg\n");
}

bool MedianFilter::test() {
    // tested by percentile filter
    return true;
}

void MedianFilter::parse(vector<string> args) {
    assert(args.size() == 1, "-medianfilter takes one argument\n");
    int radius = readInt(args[0]);
    assert(radius > -1, "radius must be positive");
    Image im = apply(stack(0), radius);
    pop();
    push(im);
}

Image MedianFilter::apply(Image im, int radius) {
    return PercentileFilter::apply(im, radius, 0.5);
}

void PercentileFilter::help() {
    printf("-percentilefilter selects a given statistical percentile over a circular support\n"
           "around each pixel. The two arguments are the support radius, and the percentile.\n"
           "A percentile argument of 0.5 gives a median filter, whereas 0 or 1 give min or\n"
           "max filters.\n\n"
           "Usage: ImageStack -load input.jpg -percentilefilter 10 0.25 -save dark.jpg\n\n");
}

bool PercentileFilter::test() {
    Image a(1024, 1024, 1, 1);
    Noise::apply(a, 0, 2);
    Image b = PercentileFilter::apply(a, 5, 0.75);
    Stats s(b);
    return nearlyEqual(s.mean(), 1.5) && nearlyEqual(s.variance(), 0);
}

void PercentileFilter::parse(vector<string> args) {
    assert(args.size() == 2, "-percentilefilter takes two arguments\n");
    int radius = readInt(args[0]);
    float percentile = readFloat(args[1]);
    assert(0 <= percentile && percentile <= 1, "percentile must be between zero and one");
    if (percentile == 1) { percentile = 0.999; }
    assert(radius > -1, "radius must be positive");
    Image im = apply(stack(0), radius, percentile);
    pop();
    push(im);
}

Image PercentileFilter::apply(Image im, int radius, float percentile) {
    struct SlidingImage {

        // We'll use a pair of heap-like data structures, with a circular
        // buffer as the leaves. The internal nodes point to the smaller
        // or greater child. Each node in the buffer belongs to at most
        // one of the two heaps at any given time.

        // Buffer to contain pixel values
        vector<float> buf;

        // The pair represents:
        // 1) Index in the circular buffer of the value at this node
        // 2) How many valid children this node has. If zero, then 1) is meaningless.
        vector<pair<int, int> > minHeap, maxHeap;

        SlidingImage(int maxKey) {
            buf.resize(maxKey);
            size_t heapSize = 1;
            while (heapSize < 2*buf.size()-1) {
                // Add a new level
                heapSize += heapSize+1;
            }
            minHeap.resize(heapSize);
            maxHeap.resize(heapSize);

            for (size_t i = 0; i < heapSize; i++) {
                minHeap[i].first = 0;
                minHeap[i].second = 0;
                maxHeap[i].first = 0;
                maxHeap[i].second = 0;
            }

            // Set the initial pointers at the leaves
            for (size_t i = 0; i < buf.size(); i++) {
                minHeap[i+buf.size()-1].first = i;
                maxHeap[i+buf.size()-1].first = i;
            }
        }

        void insert(int key, float val) {
            float p = pivot();
            buf[key] = val;
            int heapIdx = key + buf.size() - 1;
            if (isEmpty() || val < p) {
                // add to the max heap
                maxHeap[heapIdx].second = 1;
                minHeap[heapIdx].second = 0;
            } else {
                // add to the min heap
                maxHeap[heapIdx].second = 0;
                minHeap[heapIdx].second = 1;
            }
            // Fix the heaps
            updateFrom(heapIdx);
        }

        void updateFrom(int pos) {
            // walk up both heaps from the same leaf fixing pointers
            int p = pos;
            while (p) {
                // Move to the parent
                p = (p-1)/2;

                // Examine both children, and update the parent accordingly
                pair<int, int> a = minHeap[p*2+1];
                pair<int, int> b = minHeap[p*2+2];
                pair<int, int> parent;
                parent.second = a.second + b.second;
                if (a.second && b.second) {
                    parent.first = (buf[a.first] < buf[b.first]) ? a.first : b.first;
                } else if (b.second) {
                    parent.first = b.first;
                } else {
                    parent.first = a.first;
                }
                if (minHeap[p] == parent) break;
                minHeap[p] = parent;
            }

            p = pos;
            while (p) {
                p = (p-1)/2;
                pair<int, int> a = maxHeap[p*2+1];
                pair<int, int> b = maxHeap[p*2+2];
                pair<int, int> parent;
                parent.second = a.second + b.second;
                if (a.second && b.second) {
                    parent.first = (buf[a.first] > buf[b.first]) ? a.first : b.first;
                } else if (b.second) {
                    parent.first = b.first;
                } else {
                    parent.first = a.first;
                }
                if (maxHeap[p] == parent) break;
                maxHeap[p] = parent;
            }
        }

        void remove(int key) {
            int heapIdx = key+buf.size()-1;
            minHeap[heapIdx].second = 0;
            maxHeap[heapIdx].second = 0;
            updateFrom(heapIdx);
        }

        void rebalance(float p) {
            int total = maxHeap[0].second + minHeap[0].second;

            int desiredMinHeapSize = clamp(int(total * (1.0f - p)), 0, total-1);

            // Make sure there aren't too few things in the maxHeap
            while (minHeap[0].second > desiredMinHeapSize) {
                // switch the smallest thing in the minHeap into the maxHeap
                int heapIdx = minHeap[0].first + (buf.size()-1);
                minHeap[heapIdx].second = 0;
                maxHeap[heapIdx].second = 1;
                updateFrom(heapIdx);
            }

            // Make sure there aren't too many things in the maxHeap
            while (minHeap[0].second < desiredMinHeapSize) {
                // Switch the largest thing in the maxHeap into the minHeap
                int heapIdx = maxHeap[0].first + (buf.size()-1);
                minHeap[heapIdx].second = 1;
                maxHeap[heapIdx].second = 0;
                updateFrom(heapIdx);
            }
        }

        bool isEmpty() {
            return ((maxHeap[0].second + minHeap[0].second) == 0);
        }

        float pivot() {
            return buf[maxHeap[0].first];
        }

        void debug() {
            int heapSize = minHeap.size();
            printf("min heap:\n");
            for (int sz = heapSize+1; sz > 1; sz /= 2) {
                for (int i = sz/2-1; i < sz-1; i++) {
                    pair<int, int> node = minHeap[i];
                    if (node.second)
                        printf("%02d ", (int)(buf[node.first]*100));
                    else
                        printf("-- ");
                }
                printf("\n");
            }
            printf("max heap:\n");
            for (int sz = heapSize+1; sz > 1; sz /= 2) {
                for (int i = sz/2-1; i < sz-1; i++) {
                    pair<int, int> node = maxHeap[i];
                    if (node.second)
                        printf("%02d ", (int)(buf[node.first]*100));
                    else
                        printf("-- ");
                }
                printf("\n");
            }
        }

    };

    Image out(im.width, im.height, im.frames, im.channels);

    // make the filter edge profile
    int d = 2*radius+1;
    vector<int> edge(d);

    for (int i = 0; i < d; i++) {
        edge[i] = (int)(sqrtf(radius*radius - (i - radius)*(i-radius)) + 0.0001f);
    }

    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                // initialize the sliding window for this scanline
                SlidingImage window(d*d);
                for (int i = 0; i < 2*radius+1; i++) {
                    int xoff = edge[i];
                    int yoff = i - radius;

                    if (y + yoff >= im.height) { break; }
                    if (y + yoff < 0) { continue; }

                    for (int j = 0; j <= xoff; j++) {
                        if (j >= im.width) { break; }
                        float val = im(j, y+yoff, t, c);
                        window.insert(i*d + j, val);
                    }
                }

                for (int x = 0; x < im.width; x++) {
                    window.rebalance(percentile);
                    //window.debug();

                    out(x, y, t, c) = window.pivot();

                    // move the support one to the right
                    for (int i = 0; i < radius*2+1; i++) {
                        int xoff = edge[i];
                        int yoff = i - radius;

                        if (y + yoff >= im.height) { break; }
                        if (y + yoff < 0) { continue; }

                        // subtract old value
                        if (x - xoff >= 0) {
                            window.remove(i*d + (x-xoff)%d);
                        }

                        // add new value
                        if (x + xoff + 1 < im.width) {
                            float val = im(x+xoff+1, y+yoff, t, c);
                            window.insert(i*d + (x+xoff+1)%d, val);
                        }
                    }
                }
            }
        }
    }

    return out;
}



void CircularFilter::help() {
    pprintf("-circularfilter convolves the image with a uniform circular kernel. It"
            "is a good approximation to out-of-focus blur. The sole argument is the"
            "radius of the filter.\n"
            "\n"
            "Usage: ImageStack -load in.jpg -circularfilter 10 -save out.jpg\n\n");
}

bool CircularFilter::test() {
    Image impulse(21, 21, 1, 3);
    impulse(10, 10, 0, 0) = 1;
    impulse(10, 10, 0, 1) = 2;
    impulse(10, 10, 0, 2) = 3;
    Image blurry = CircularFilter::apply(impulse, 5);
    float ratio = blurry(10, 10, 0, 0);
    for (int y = 0; y < 21; y++) {
        float fy = (y - 10.0f)/5;
        for (int x = 0; x < 21; x++) {
            float fx = (x - 10.0f)/5;
            float correct = (fx*fx + fy*fy) <= 1.00001 ? ratio : 0;
            if (!nearlyEqual(blurry(x, y, 0, 0), correct*1)) return false;
            if (!nearlyEqual(blurry(x, y, 0, 1), correct*2)) return false;
            if (!nearlyEqual(blurry(x, y, 0, 2), correct*3)) return false;
        }
    }
    return nearlyEqual(Stats(blurry).sum(), 6);
}

void CircularFilter::parse(vector<string> args) {
    assert(args.size() == 1, "-circularfilter takes one argument\n");

    Image im = apply(stack(0), readInt(args[0]));
    pop();
    push(im);
}

Image CircularFilter::apply(Image im, int radius) {
    Image out(im.width, im.height, im.frames, im.channels);

    // maintain the average response currently under the filter, and the number of pixels under the filter
    float average = 0;
    int count = 0;

    // make the filter edge profile
    vector<int> edge(radius*2+1);
    for (int i = 0; i < 2*radius+1; i++) {
        edge[i] = (int)(sqrtf(radius*radius - (i - radius)*(i-radius)) + 0.0001f);
    }

    // figure out the filter area
    for (int i = 0; i < 2*radius+1; i++) {
        count += edge[i]*2+1;
    }

    float invArea = 1.0f/count;

    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                average = 0;
                // initialize the average and count
                for (int i = 0; i < 2*radius+1; i++) {
                    int xoff = edge[i];
                    int yoff = i - radius;
                    int realY = clamp(y + yoff, 0, im.height-1);

                    for (int x = -xoff; x <= xoff; x++) {
                        int realX = clamp(x, 0, im.width-1);
                        float val = im(realX, realY, t, c);
                        average += val;
                    }
                }

                for (int x = 0; x < im.width; x++) {
                    out(x, y, t, c) = average * invArea;

                    // move the histogram to the right
                    for (int i = 0; i < radius*2+1; i++) {
                        int realXOld = max(0, x-edge[i]);
                        int realXNew = min(x+edge[i]+1, im.width-1);
                        int realY = clamp(y+i-radius, 0, im.height-1);

                        // add new value, subtract old value
                        average += im(realXNew, realY, t, c);
                        average -= im(realXOld, realY, t, c);
                    }
                }
            }
        }
    }

    return out;
}



void Envelope::help() {
    pprintf("-envelope computes a lower or upper envelope of the input, which is"
            " smooth, and less than (or greater than) the input. The first argument"
            " should be \"lower\" or \"upper\". The second argument is the desired"
            " smoothness, which is roughly proportional to the pixel radius of a blur.\n"
            "\n"
            "Usage: ImageStack -load a.jpg -envelope upper 50 -display\n");
}

bool Envelope::test() {
    {
        Image a(123, 234, 2, 3);
        Noise::apply(a, 0, 1);
        Image b = a.copy();
        Envelope::apply(b, Upper, 3);
        b -= a;
        if (!nearlyEqual(Stats(b).minimum(), 0)) return false;
    }

    {
        Image a(123, 234, 2, 3);
        Noise::apply(a, 0, 1);
        Image b = a.copy();
        Envelope::apply(b, Lower, 3);
        a -= b;
        if (!nearlyEqual(Stats(a).minimum(), 0)) return false;
    }

    return true;
}

void Envelope::parse(vector<string> args) {
    assert(args.size() == 2, "-envelope takes two arguments\n");
    Mode m = Lower;
    if (args[0] == "lower") { m = Lower; }
    else if (args[0] == "upper") { m = Upper; }
    else { panic("Unknown mode: %s. Must be lower or upper.\n", args[0].c_str()); }

    apply(stack(0), m, readInt(args[1]));
}

void Envelope::apply(Image im, Mode m, int radius) {
    if (m == Upper) {
        MaxFilter::apply(im, radius);
        RectFilter::apply(im, 2*radius+1, 2*radius+1, 1);
        radius = (radius+2)/3;
        MaxFilter::apply(im, radius);
        RectFilter::apply(im, 2*radius+1, 2*radius+1, 1);
    }

    if (m == Lower) {
        MinFilter::apply(im, radius);
        RectFilter::apply(im, 2*radius+1, 2*radius+1, 1);
        radius = (radius+2)/3;
        MinFilter::apply(im, radius);
        RectFilter::apply(im, 2*radius+1, 2*radius+1, 1);
    }
}

void HotPixelSuppression::help() {
    pprintf("-hotpixelsuppression removes salt-and-pepper noise from an image by"
            " constraining each pixel to be within the bounds of its four"
            " neighbors\n\n"
            "Usage: ImageStack -load noisy.jpg -hotpixelsuppression -save denoised.jpg\n");
}

bool HotPixelSuppression::test() {
    Image a(100, 100, 1, 1);
    Noise::apply(a, -4, 12);
    Image b = HotPixelSuppression::apply(a);
    for (int y = 1; y < 99; y++) {
        for (int x = 1; x < 99; x++) {
            if (b(x, y) > a(x-1, y) &&
                b(x, y) > a(x+1, y) &&
                b(x, y) > a(x, y-1) &&
                b(x, y) > a(x, y+1)) return false;
            if (b(x, y) < a(x-1, y) &&
                b(x, y) < a(x+1, y) &&
                b(x, y) < a(x, y-1) &&
                b(x, y) < a(x, y+1)) return false;
        }
    }
    return true;
}

void HotPixelSuppression::parse(vector<string> args) {
    assert(args.size() == 0,
           "-hotpixelsuppression takes no arguments\n");
    Image im = apply(stack(0));
    pop();
    push(im);
}

Image HotPixelSuppression::apply(Image im) {
    Image out(im.width, im.height, im.frames, im.channels);

    for (int t = 0; t < im.frames; t++) {
        for (int y = 0; y < im.height; y++) {
            for (int x = 0; x < im.width; x++) {
                for (int c = 0; c < im.channels; c++) {
                    float maxn = -INF;
                    float minn = INF;
                    if (x > 0) {
                        maxn = max(maxn, im(x-1, y, t, c));
                        minn = min(minn, im(x-1, y, t, c));
                    }
                    if (x < im.width-1) {
                        maxn = max(maxn, im(x+1, y, t, c));
                        minn = min(minn, im(x+1, y, t, c));
                    }
                    if (y > 0) {
                        maxn = max(maxn, im(x, y-1, t, c));
                        minn = min(minn, im(x, y-1, t, c));
                    }
                    if (y < im.height-1) {
                        maxn = max(maxn, im(x, y+1, t, c));
                        minn = min(minn, im(x, y+1, t, c));
                    }
                    float here = im(x, y, t, c);
                    if (here > maxn) here = maxn;
                    if (here < minn) here = minn;
                    out(x, y, t, c) = here;
                }
            }
        }
    }

    return out;
}

#include "footer.h"
