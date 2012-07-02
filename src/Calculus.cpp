#include "main.h"
#include "Calculus.h"
#include "Geometry.h"
#include "Arithmetic.h"
#include "Stack.h"
#include "Convolve.h"
#include "Filter.h"
#include "File.h"
#include "Display.h"
#include "LAHBPCG.h"
#include "Statistics.h"
#include "header.h"

void Gradient::help() {
    pprintf("-gradient takes the backward differences in the dimension specified by"
            " the argument. Values outside the image are assumed to be zero, so the"
            " first row, or column, or frame, will not change, effectively storing"
            " the initial value to make later integration easy. Multiple arguments"
            " can be given to differentiate with respect to multiple dimensions in"
            " order (although the order does not matter).\n"
            "\n"
            "Warning: Don't expect to differentiate more than twice and be able to"
            " get back the image by integrating. Numerical errors will dominate.\n"
            "\n"
            "Usage: ImageStack -load a.jpg -gradient x y -save out.jpg\n");
}

bool Gradient::test() {
    Image a(16, 14, 18, 3);
    Noise::apply(a, 0, 1);
    Image dx = a.copy();
    Gradient::apply(dx, 'x');
    Image dy = a.copy();
    Gradient::apply(dy, 'y');
    Image dt = a.copy();
    Gradient::apply(dt, 't');
    for (int i = 0; i < 100; i++) {
        int x = randomInt(1, a.width-1);
        int y = randomInt(1, a.height-1);
        int t = randomInt(1, a.frames-1);
        int c = randomInt(1, a.channels-1);
        if (!nearlyEqual(dx(x, y, t, c), a(x, y, t, c) - a(x-1, y, t, c))) return false;
        if (!nearlyEqual(dy(x, y, t, c), a(x, y, t, c) - a(x, y-1, t, c))) return false;
        if (!nearlyEqual(dt(x, y, t, c), a(x, y, t, c) - a(x, y, t-1, c))) return false;
    }
    return true;
}

void Gradient::parse(vector<string> args) {
    assert(args.size() > 0, "-gradient requires at least one argument\n");
    for (size_t i = 0; i < args.size(); i++) {
        apply(stack(0), args[i]);
    }
}

// gradient can be called as gradient('t') or gradient("xyt")
void Gradient::apply(Image im, string dimensions) {
    for (size_t i = 0; i < dimensions.size(); i++) {
        apply(im, dimensions[i]);
    }
}

void Gradient::apply(Image im, char dimension) {
    int mint = 0, minx = 0, miny = 0;
    int dt = 0, dx = 0, dy = 0;

    if (dimension == 'x') {
        dx = 1;
        minx = 1;
    } else if (dimension == 'y') {
        dy = 1;
        miny = 1;
    } else if (dimension == 't') {
        dt = 1;
        mint = 1;
    } else {
        panic("Must differentiate with respect to x, y, or t\n");
    }

    // walk backwards through the data, looking at the untouched data for the differences
    for (int c = 0; c < im.channels; c++) {
        for (int t = im.frames - 1; t >= mint; t--) {
            for (int y = im.height - 1; y >= miny; y--) {
                for (int x = im.width - 1; x >= minx; x--) {
                    im(x, y, t, c) -= im(x - dx, y - dy, t - dt, c);
                }
            }
        }
    }
}


void Integrate::help() {
    pprintf("-integrate computes partial sums along the given dimension. It is the"
            " inverse of the -gradient operator. Multiply dimensions can be given"
            " as arguments, for example -integrate x y will produce a summed area"
            " table of an image. Allowed dimensions are x, y, or t.\n"
            "\n"
            "Warning: Don't expect to integrate more than twice and be able to get"
            " back the image by differentiating. Numerical errors will dominate.\n"
            "\n"
            "Usage: ImageStack -load a.jpg -gradient x y -integrate y x -save a.jpg\n");
}

bool Integrate::test() {
    Image a(23, 13, 32, 3);
    Noise::apply(a, 0, 1);
    Image b = a.copy();
    Gradient::apply(b, "xyt");
    Integrate::apply(b, "xty");
    return nearlyEqual(a, b);
}

void Integrate::parse(vector<string> args) {
    assert(args.size() > 0, "-integrate requires at least one argument\n");
    for (size_t i = 0; i < args.size(); i++) {
        apply(stack(0), args[i]);
    }
}

// integrate can be called as integrate('t') or integrate("xyt")
void Integrate::apply(Image im, string dimensions) {
    for (size_t i = 0; i < dimensions.size(); i++) {
        apply(im, dimensions[i]);
    }
}

void Integrate::apply(Image im, char dimension) {
    int minx = 0, miny = 0, mint = 0;
    int dx = 0, dy = 0, dt = 0;

    if (dimension == 'x') {
        dx = 1;
        minx = 1;
    } else if (dimension == 'y') {
        dy = 1;
        miny = 1;
    } else if (dimension == 't') {
        dt = 1;
        mint = 1;
    } else {
        panic("Must integrate with respect to x, y, or t\n");
    }

    // walk forwards through the data, adding up as we go
    for (int c = 0; c < im.channels; c++) {
        for (int t = mint; t < im.frames; t++) {
            for (int y = miny; y < im.height; y++) {
                for (int x = minx; x < im.width; x++) {
                    im(x, y, t, c) += im(x - dx, y - dy, t - dt, c);
                }
            }
        }
    }

}


void GradMag::help() {
    pprintf("-gradmag computes the square gradient magnitude at each pixel in x and"
            " y. Temporal gradients are ignored. The gradient is estimated using"
            " backward differences, and the image is assumed to be zero outside its"
            " bounds.\n"
            "\n"
            "Usage: ImageStack -load input.jpg -gradmag -save out.jpg\n");
}

bool GradMag::test() {
    Image a(233, 123, 1, 5);
    Noise::apply(a, 0, 1);
    Image dx = a.copy();
    Gradient::apply(dx, 'x');
    Image dy = a.copy();
    Gradient::apply(dy, 'y');
    GradMag::apply(a);
    return nearlyEqual(a, dx*dx + dy*dy);
}


void GradMag::parse(vector<string> args) {
    assert(args.size() == 0, "-gradmag takes no arguments\n");
    apply(stack(0));
}

void GradMag::apply(Image im) {
    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames; t++) {
            for (int y = im.height-1; y >=0; y--) {
                for (int x = im.width-1; x >= 0; x--) {
                    float dx = im(x, y, t, c) - (x > 0 ? im(x-1, y, t, c) : 0);
                    float dy = im(x, y, t, c) - (y > 0 ? im(x, y-1, t, c) : 0);
                    im(x, y, t, c) = dx*dx + dy*dy;
                }
            }
        }
    }
}


void Poisson::help() {
    pprintf("-poisson assumes the stack contains gradients images in x and y, and"
            " attempts to find the image which fits those gradients best in a least"
            " squares sense. It uses a preconditioned conjugate gradient descent"
            " method. It takes one argument, which is required RMS error of the"
            " result. This defaults to 0.01 if not given.\n"
            "\n"
            "Usage: ImageStack -load dx.tmp -load dy.tmp \n"
            "                  -poisson 0.0001 -save out.jpg\n\n");
}

bool Poisson::test() {
    Image a(233, 123, 1, 5);
    Noise::apply(a, 0, 1);
    Image dx = a.copy();
    Gradient::apply(dx, 'x');
    Image dy = a.copy();
    Gradient::apply(dy, 'y');
    Image b = Poisson::apply(dx, dy, 0.00001);
    return nearlyEqual(a, b);
}

void Poisson::parse(vector<string> args) {
    assert(args.size() < 2, "-poisson requires one or fewer arguments\n");
    float rms = 0.01;
    if (args.size() > 0) {
        rms = readFloat(args[0]);
    }

    push(apply(stack(1), stack(0), rms));
}

Image Poisson::apply(Image dx, Image dy, float rms) {
    assert(dx.width  == dy.width &&
           dx.height == dy.height &&
           dx.frames == dy.frames &&
           dx.channels == dy.channels,
           "derivatives must be matching size and number of channels\n");


    Image zerosc(dx.width, dx.height, dx.frames, dx.channels);
    Image zeros1(dx.width, dx.height, dx.frames, 1);
    Image ones1(dx.width, dx.height, dx.frames, 1);
    ones1.set(1);
    return LAHBPCG::apply(zerosc, dx, dy, zeros1, ones1, ones1, 999999, rms);
}

#include "footer.h"
