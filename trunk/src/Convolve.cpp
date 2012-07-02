#include "main.h"
#include "Convolve.h"
#include "Arithmetic.h"
#include "Geometry.h"
#include "DFT.h"
#include "File.h"
#include "Statistics.h"
#include "header.h"

void Convolve::help() {
    pprintf("-convolve takes a width, height, and frames and a single-channel 3D"
            " kernel specified across the rows, then down the columns, then over"
            " time, and convolves the current image by that matrix independently in"
            " each channel.\n"
            "\n"
            "With no numeric arguments, -convolve will use the next image on the stack as"
            " the filter.\n"
            "\n"
            "Boundary conditions can be specified by appending the argument"
            " \"zero\", \"homogeneous\", \"clamp\", or \"wrap\", which result in the"
            " respective assumptions: the image is zero outside the boundary; the"
            " image is weighted with weight one inside the boundary, and weight zero"
            " outside the boundary; the image values outside the boundary clamp to"
            " the nearest defined image value (Neumann); and the image wraps outside"
            " the boundary.\n"
            "\n"
            "Convolution by multi-channel filters is poorly defined, because it"
            " requires a vector-vector multiplication between filter values and"
            " image values. By specifying a final argument of \"inner\", \"outer\","
            " or \"elementwise\", the multiplication used is correspondingly the"
            " inner product (or matrix product if the image and kernel have a"
            " differing number of channels); the outer product; or an elementwise"
            " product. If the kernel has k channels and the image has m channels,"
            " \"inner\" produces an image with max(m/k, k/m) channels, \"outer\""
            " produces an image with m*k channels, and \"elementwise\" requires"
            " that m==k and produces an image with the same number of channels. The"
            " default method is \"outer\" if the channel counts are different, and "
            "\"elementwise\" if they are the same.\n"
            "\n"
            "Taking a horizontal gradient with zero boundary condition: \n"
            " ImageStack -load a.tga -convolve 2 1 1  -1 1 zero -save dx.tga\n"
            "Convolving by a bank of filters: \n"
            " ImageStack -load bank.tmp -load a.tga -convolve homogeneous outer\n");
}

bool Convolve::test() {
    Image impulse(32, 32, 32, 2);
    impulse(15, 15, 15, 0) = 1;
    impulse(15, 15, 15, 1) = 2;
    Image kernel(5, 5, 5, 4);
    Noise::apply(kernel, 0, 10);
    Image correct(32, 32, 32, 2);
    correct
    .region(13, 13, 13, 0, 5, 5, 5, 1)
    .set(1*kernel.channel(0) + 2*kernel.channel(1));
    correct
    .region(13, 13, 13, 1, 5, 5, 5, 1)
    .set(1*kernel.channel(2) + 2*kernel.channel(3));
    Image result = Convolve::apply(impulse, kernel, Zero, Multiply::Inner);
    return nearlyEqual(result, correct);
}

void Convolve::parse(vector<string> args) {
    string boundaryCondition = "homogeneous";
    string channelMode = "outer";

    Image filter;

    if (args.size() > 3) {
        int frames, width, height;
        size_t size;
        width = readInt(args[0]);
        height = readInt(args[1]);
        frames = readInt(args[2]);
        size = frames * width * height;
        assert(args.size() >= size+3,
               "a size of %ix%ix%i requires at least %i more arguments. %i were given.",
               width, height, frames, size, (int)args.size() - 3);
        assert(size % 2 == 1, "filter must have odd size\n");

        assert(args.size() <= size+4,
               "a size of %ix%ix%i requires at most %i more arguments. %i were given.",
               width, height, frames, size, (int)args.size() - 2);

        filter = Image(width, height, frames, 1);

        size_t i = 3;
        for (int t = 0; t < frames; t++) {
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    filter(x, y, t, 0) = readFloat(args[i++]);
                }
            }
        }

        if (args.size() == size+4) {
            boundaryCondition = args[size+3];
        }

    } else if (args.size() < 3) {
        filter = stack(1);
        if (args.size() >= 1) {
            boundaryCondition = args[0];
        }
        if (args.size() == 2) {
            channelMode = args[1];
        } else if (stack(0).channels == filter.channels) {
            channelMode = "elementwise";
        }
    } else {
        panic("-convolve needs either zero, one, two, or at least four arguments\n");
    }

    Multiply::Mode m = Multiply::Outer;
    BoundaryCondition b = Homogeneous;

    if (boundaryCondition == "zero") { b = Zero; }
    else if (boundaryCondition == "homogeneous") { b = Homogeneous; }
    else if (boundaryCondition == "clamp") { b = Clamp; }
    else if (boundaryCondition == "wrap") { b = Wrap; }
    else {
        panic("Unknown boundary condition: %s\n", boundaryCondition.c_str());
    }

    if (channelMode == "inner") { m = Multiply::Inner; }
    else if (channelMode == "outer") { m = Multiply::Outer; }
    else if (channelMode == "elementwise") { m = Multiply::Elementwise; }
    else {
        panic("Unknown vector-vector multiplication: %s\n", channelMode.c_str());
    }

    Image im = apply(stack(0), filter, b, m);
    pop();
    push(im);

}

// For a single channel, out += in * filter
void Convolve::convolveSingle(Image in, Image filter, Image out,
                              BoundaryCondition b) {
    assert(in.channels == 1 && filter.channels == 1 && out.channels == 1,
           "convolveSingle should only be called on single-channel images");

    int filterSize = filter.frames * filter.width * filter.height;
    assert(filterSize % 2 == 1, "filter must have odd size (%d %d %d)\n", filter.width, filter.height, filter.frames);

    int xoff = (filter.width - 1)/2;
    int yoff = (filter.height - 1)/2;
    int toff = (filter.frames - 1)/2;

    if (b == Zero) {
        for (int t = 0; t < in.frames; t++) {
            for (int y = 0; y < in.height; y++) {
                for (int x = 0; x < in.width; x++) {
                    float v = 0;
                    for (int dt = -toff; dt <= toff; dt++) {
                        if (t + dt < 0) continue;
                        if (t + dt >= in.frames) break;
                        for (int dy = -yoff; dy <= yoff; dy++) {
                            if (y + dy < 0) continue;
                            if (y + dy >= in.height) break;
                            for (int dx = -xoff; dx <= xoff; dx++) {
                                if (x + dx < 0) continue;
                                if (x + dx >= in.width) break;
                                float w = filter(xoff-dx, yoff-dy, toff-dt, 0);
                                v += in(x+dx, y+dy, t+dt, 0) * w;
                            }
                        }
                    }
                    out(x, y, t, 0) += v;
                }
            }
        }
    } else if (b == Homogeneous) {
        float filterSum = Stats(filter).sum();
        for (int t = 0; t < in.frames; t++) {
            for (int y = 0; y < in.height; y++) {
                for (int x = 0; x < in.width; x++) {
                    float weightSum = 0;
                    float v = 0;
                    for (int dt = -toff; dt <= toff; dt++) {
                        if (t + dt < 0) continue;
                        if (t + dt >= in.frames) break;
                        for (int dy = -yoff; dy <= yoff; dy++) {
                            if (y + dy < 0) continue;
                            if (y + dy >= in.height) break;
                            for (int dx = -xoff; dx <= xoff; dx++) {
                                if (x + dx < 0) continue;
                                if (x + dx >= in.width) break;
                                float w = filter(xoff-dx, yoff-dy, toff-dt, 0);
                                v += in(x+dx, y+dy, t+dt, 0) * w;
                                weightSum += w;
                            }
                        }
                    }
                    if (filterSum != weightSum) {
                        v *= filterSum / weightSum;
                    }
                    out(x, y, t, 0) += v;
                }
            }
        }
    } else if (b == Clamp) {
        for (int t = 0; t < in.frames; t++) {
            for (int y = 0; y < in.height; y++) {
                for (int x = 0; x < in.width; x++) {
                    float v = 0;
                    for (int dt = -toff; dt <= toff; dt++) {
                        int tc = clamp(t+dt, 0, in.frames-1);
                        for (int dy = -yoff; dy <= yoff; dy++) {
                            int yc = clamp(y+dy, 0, in.height-1);
                            for (int dx = -xoff; dx <= xoff; dx++) {
                                int xc = clamp(x+dx, 0, in.width-1);
                                float w = filter(xoff-dx, yoff-dy, toff-dt, 0);
                                v += in(xc, yc, tc, 0) * w;
                            }
                        }
                    }
                    out(x, y, t, 0) += v;
                }
            }
        }
    } else if (b == Wrap) {
        for (int t = 0; t < in.frames; t++) {
            for (int y = 0; y < in.height; y++) {
                for (int x = 0; x < in.width; x++) {
                    float v = 0;
                    for (int dt = -toff; dt <= toff; dt++) {
                        int tc = (t+dt+toff*in.frames)%in.frames;
                        for (int dy = -yoff; dy <= yoff; dy++) {
                            int yc = (y+dy+yoff*in.height)%in.height;
                            for (int dx = -xoff; dx <= xoff; dx++) {
                                int xc = (x+dx+xoff*in.width)%in.width;
                                float w = filter(xoff-dx, yoff-dy, toff-dt, 0);
                                v += in(xc, yc, tc, 0) * w;
                            }
                        }
                    }
                    out(x, y, t, 0) += v;
                }
            }
        }
    } else {
        panic("Unknown boundary condition");
    }
}

Image Convolve::apply(Image im, Image filter, BoundaryCondition b, Multiply::Mode m) {
    Image out;
    if (m == Multiply::Inner) {
        assert(filter.channels % im.channels == 0 ||
               im.channels % filter.channels == 0,
               "To perform an inner or matrix product, the channel count "
               "of either the image or the filter must be a multiple of "
               "the channel count of the other.");
        if (im.channels < filter.channels) {
            out = Image(im.width, im.height, im.frames, filter.channels/im.channels);
            for (int i = 0; i < filter.channels; i++) {
                convolveSingle(im.channel(i % im.channels),
                               filter.channel(i),
                               out.channel(i / im.channels), b);
            }
        } else {
            out = Image(im.width, im.height, im.frames, im.channels/filter.channels);
            for (int i = 0; i < im.channels; i++) {
                convolveSingle(im.channel(i),
                               filter.channel(i % filter.channels),
                               out.channel(i / filter.channels), b);
            }
        }
    } else if (m == Multiply::Outer) {
        out = Image(im.width, im.height, im.frames, im.channels * filter.channels);
        for (int i = 0; i < im.channels; i++) {
            for (int j = 0; j < filter.channels; j++) {
                convolveSingle(im.channel(i),
                               filter.channel(j),
                               out.channel(i*filter.channels + j), b);
            }
        }
    } else if (m == Multiply::Elementwise) {
        assert(im.channels == filter.channels,
               "For element-wise multiplication, the image "
               "and filter must have the same number of channels.");
        out = Image(im.width, im.height, im.frames, im.channels);
        for (int i = 0; i < im.channels; i++) {
            convolveSingle(im.channel(i), filter.channel(i), out.channel(i), b);
        }

    } else {
        panic("Unknown multiplication mode");
    }
    return out;
}

#include "footer.h"

