#include "main.h"
#include "Wavelet.h"
#include "Geometry.h"
#include "Statistics.h"
#include "Calculus.h"
#include "Convolve.h"
#include "header.h"

void Haar::help() {
    pprintf("-haar performs the standard 2D haar transform of an image. The image"
            " size must be a power of two. If given an integer argument k, it only"
            " recurses k times, and the image size must be a multiple of 2^k.\n"
            "\n"
            "Usage: ImageStack -load in.jpg -haar 1 -save out.jpg\n\n");
}

bool Haar::test() {
    Image a(256, 256, 3, 1);
    Noise::apply(a, -34, 23);
    Image b = a.copy();
    Haar::apply(b, 4);

    // Top left 16x16 should be a downsampled version of the original
    Image smaller = Downsample::apply(a, 16, 16, 1);
    if (!nearlyEqual(smaller, b.region(0, 0, 0, 0, 16, 16, 3, 1))) return false;

    // Bottom right 128x128 should be a subsampled derivative
    Image deriv = a.copy();
    Gradient::apply(deriv, "xy");
    deriv = Subsample::apply(deriv, 2, 2, 1, 1);
    return nearlyEqual(deriv, b.region(128, 128, 0, 0, 128, 128, 3, 1));
}

void Haar::parse(vector<string> args) {
    if (args.size() == 0) {
        apply(stack(0));
    } else if (args.size() == 1) {
        apply(stack(0), readInt(args[0]));
    } else {
        panic("-haar requires zero or one arguments\n");
    }
}

void Haar::apply(Image im, int times) {

    if (times <= 0) {
        assert(im.width == im.height, "to perform a full haar transform, the image must be square.\n");
        times = 0;
        int w = im.width >> 1;
        while (w) {
            times++;
            w >>= 1;
        }
    }

    int factor = 1 << times;
    assert(im.width % factor == 0, "the image width is not a multiple of 2^%i", times);
    assert(im.height % factor == 0, "the image height is not a multiple of 2^%i", times);

    // transform in x
    Image win = im.region(0, 0, 0, 0, im.width, im.height, im.frames, im.channels);
    for (int i = 0; i < times; i++) {
        for (int c = 0; c < win.channels; c++) {
            for (int t = 0; t < win.frames; t++) {
                for (int y = 0; y < win.height; y++) {
                    for (int x = 0; x < win.width; x+=2) {
                        float a = win(x, y, t, c);
                        float b = win(x+1, y, t, c);
                        win(x, y, t, c) = (a+b)/2;
                        win(x+1, y, t, c) = b-a;
                    }
                }
            }
        }
        // separate into averages and differences
        Deinterleave::apply(win, 2, 1, 1);
        // repeat on the averages
        win = win.region(0, 0, 0, 0, win.width/2, win.height, win.frames, win.channels);
    }

    // transform in y
    win = im.region(0, 0, 0, 0, im.width, im.height, im.frames, im.channels);
    for (int i = 0; i < times; i++) {
        for (int c = 0; c < win.channels; c++) {
            for (int t = 0; t < win.frames; t++) {
                for (int y = 0; y < win.height; y+=2) {
                    for (int x = 0; x < win.width; x++) {
                        float a = win(x, y, t, c);
                        float b = win(x, y+1, t, c);
                        win(x, y, t, c) = (a+b)/2;
                        win(x, y+1, t, c) = b-a;
                    }
                }
            }
        }
        // separate into averages and differences
        Deinterleave::apply(win, 1, 2, 1);
        // repeat on the averages
        win = win.region(0, 0, 0, 0, win.width, win.height/2, win.frames, win.channels);
    }
}


void InverseHaar::help() {
    pprintf("-inversehaar inverts the haar transformation with the same"
            " argument. See -help haar for detail.\n");

}

bool InverseHaar::test() {
    Image a(256, 256, 3, 3);
    Noise::apply(a, 0, 1);
    Image b = a.copy();
    Haar::apply(b, 3);
    InverseHaar::apply(b, 3);
    return nearlyEqual(a, b);
}

void InverseHaar::parse(vector<string> args) {
    if (args.size() == 0) {
        apply(stack(0));
    } else if (args.size() == 1) {
        apply(stack(0), readInt(args[0]));
    } else {
        panic("-haar requires zero or one arguments\n");
    }
}

void InverseHaar::apply(Image im, int times) {
    if (times <= 0) {
        assert(im.width == im.height, "to perform a full haar transorm, the image must be square.\n");
        times = 0;
        int w = im.width >> 1;
        while (w) {
            times++;
            w >>= 1;
        }
    }

    int factor = 1 << times;
    assert(im.width % factor == 0, "the image width is not a multiple of 2^%i", times);
    assert(im.height % factor == 0, "the image height is not a multiple of 2^%i", times);

    // transform in y
    int h = 2*im.height/factor;
    Image win = im.region(0, 0, 0, 0, im.width, h, im.frames, im.channels);
    while (1) {
        // combine the averages and differences
        Interleave::apply(win, 1, 2, 1);

        for (int c = 0; c < win.frames; c++) {
            for (int t = 0; t < win.frames; t++) {
                for (int y = 0; y < win.height; y+=2) {
                    for (int x = 0; x < win.width; x++) {
                        float avg = win(x, y, t, c);
                        float diff = win(x, y+1, t, c);
                        win(x, y, t, c) = avg-diff/2;
                        win(x, y+1, t, c) = avg+diff/2;
                    }
                }
            }
        }
        // repeat
        h *= 2;
        if (h > im.height) { break; }
        win = im.region(0, 0, 0, 0, im.width, h, im.frames, im.channels);
    }

    // transform in x
    int w = 2*im.width/factor;
    win = im.region(0, 0, 0, 0, w, im.height, im.frames, im.channels);
    while (1) {
        // combine the averages and differences
        Interleave::apply(win, 2, 1, 1);

        for (int c = 0; c < win.channels; c++) {
            for (int t = 0; t < win.frames; t++) {
                for (int y = 0; y < win.height; y++) {
                    for (int x = 0; x < win.width; x+=2) {
                        float avg = win(x, y, t, c);
                        float diff = win(x+1, y, t, c);
                        win(x, y, t, c) = avg-diff/2;
                        win(x+1, y, t, c) = avg+diff/2;
                    }
                }
            }
        }
        // repeat
        w *= 2;
        if (w > im.width) { break; }
        win = im.region(0, 0, 0, 0, w, im.height, im.frames, im.channels);
    }
}




#define DAUB0 0.4829629131445341
#define DAUB1 0.83651630373780772
#define DAUB2 0.22414386804201339
#define DAUB3 -0.12940952255126034

void Daubechies::help() {
    pprintf("-daubechies performs the standard 2D daubechies 4 wavelet transform of"
            " an image. The image size must be a power of two.\n"
            "\n"
            "Usage: ImageStack -load in.jpg -daubechies -save out.jpg\n\n");
}

bool Daubechies::test() {
    Image a(256, 256, 3, 3);
    Noise::apply(a, 0, 1);
    Image filter(7, 1, 1, 1);
    filter(0, 0) = -DAUB0;
    filter(1, 0) = DAUB1;
    filter(2, 0) = -DAUB2;
    filter(3, 0) = DAUB3;
    filter(4, 0) = 0;
    filter(5, 0) = 0;
    filter(6, 0) = 0;
    Image b = Convolve::apply(a, filter);
    b = Convolve::apply(b, Transpose::apply(filter, 'x', 'y'));
    b = Subsample::apply(b, 2, 2, 0, 0);
    Daubechies::apply(a);
    return nearlyEqual(b.region(0, 0, 0, 0, 100, 100, 3, 3),
                       a.region(128, 128, 0, 0, 100, 100, 3, 3));
}

void Daubechies::parse(vector<string> args) {
    assert(args.size() == 0, "-daubechies takes no arguments");
    apply(stack(0));
}

void Daubechies::apply(Image im) {

    int i;
    for (i = 1; i < im.width; i <<= 1);
    assert(i == im.width, "Image width must be a power of two\n");
    for (i = 1; i < im.height; i <<= 1);
    assert(i == im.height, "Image height must be a power of two\n");

    // transform in x
    for (int width = im.width; width > 1; width /= 2) {
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int y = 0; y < im.height; y++) {
                    float saved1st = im(0, y, t, c);
                    float saved2nd = im(1, y, t, c);

                    for (int x = 0; x < width-2; x+=2) {
                        float v0 = im(x, y, t, c);
                        float v1 = im(x+1, y, t, c);
                        float v2 = im(x+2, y, t, c);
                        float v3 = im(x+3, y, t, c);
                        im(x, y, t, c) = DAUB0 * v0 + DAUB1 * v1 + DAUB2 * v2 + DAUB3 * v3;
                        im(x+1, y, t, c) = DAUB3 * v0 - DAUB2 * v1 + DAUB1 * v2 - DAUB0 * v3;
                    }
                    // special case the last two elements using rotation
                    float v0 = im(width-2, y, t, c);
                    float v1 = im(width-1, y, t, c);
                    float v2 = saved1st;
                    float v3 = saved2nd;
                    im(width-2, y, t, c) = DAUB0 * v0 + DAUB1 * v1 + DAUB2 * v2 + DAUB3 * v3;
                    im(width-1, y, t, c) = DAUB3 * v0 - DAUB2 * v1 + DAUB1 * v2 - DAUB0 * v3;
                }
            }
        }
        // separate into averages and differences
        Deinterleave::apply(im.region(0, 0, 0, 0,
                                      width, im.height, im.frames, im.channels),
                            2, 1, 1);
    }

    // transform in y
    for (int height = im.height; height > 1; height /= 2) {
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int x = 0; x < im.width; x++) {
                    float saved1st = im(x, 0, t, c);
                    float saved2nd = im(x, 1, t, c);
                    for (int y = 0; y < height-2; y+=2) {
                        float v0 = im(x, y, t, c);
                        float v1 = im(x, y+1, t, c);
                        float v2 = im(x, y+2, t, c);
                        float v3 = im(x, y+3, t, c);
                        im(x, y, t, c) = DAUB0 * v0 + DAUB1 * v1 + DAUB2 * v2 + DAUB3 * v3;
                        im(x, y+1, t, c) = DAUB3 * v0 - DAUB2 * v1 + DAUB1 * v2 - DAUB0 * v3;
                    }
                    // special case the last two elements using rotation
                    float v0 = im(x, height-2, t, c);
                    float v1 = im(x, height-1, t, c);
                    float v2 = saved1st;
                    float v3 = saved2nd;
                    im(x, height-2, t, c) = DAUB0 * v0 + DAUB1 * v1 + DAUB2 * v2 + DAUB3 * v3;
                    im(x, height-1, t, c) = DAUB3 * v0 - DAUB2 * v1 + DAUB1 * v2 - DAUB0 * v3;
                }
            }
        }
        // separate into averages and differences
        Deinterleave::apply(im.region(0, 0, 0, 0,
                                      im.width, height, im.frames, im.channels),
                            1, 2, 1);
    }
}


void InverseDaubechies::help() {
    pprintf("-inversedaubechies performs the standard 2D daubechies 4 wavelet "
            " transform of an image. The image size must be a power of two.\n"
            "\n"
            "Usage: ImageStack -load in.jpg -inversedaubechies -save out.jpg\n");
}

bool InverseDaubechies::test() {
    Image a(256, 256, 3, 3);
    Noise::apply(a, -5, 12);
    Image b = a.copy();
    Daubechies::apply(b);
    InverseDaubechies::apply(b);
    return nearlyEqual(a, b);
}

void InverseDaubechies::parse(vector<string> args) {
    assert(args.size() == 0, "-inversedaubechies takes no arguments");
    apply(stack(0));
}

void InverseDaubechies::apply(Image im) {

    int i;
    for (i = 1; i < im.width; i <<= 1) { ; }
    assert(i == im.width, "Image width must be a power of two\n");
    for (i = 1; i < im.height; i <<= 1) { ; }
    assert(i == im.height, "Image height must be a power of two\n");


    // transform in x
    for (int width = 2; width <= im.width; width *= 2) {
        // Collect averages and differences
        Interleave::apply(im.region(0, 0, 0, 0,
                                    width, im.height, im.frames, im.channels),
                          2, 1, 1);

        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int y = 0; y < im.height; y++) {
                    float saved1st = im(width-1, y, t, c);
                    float saved2nd = im(width-2, y, t, c);

                    for (int x = width-4; x >= 0; x-=2) {
                        float v0 = im(x, y, t, c);
                        float v1 = im(x+1, y, t, c);
                        float v2 = im(x+2, y, t, c);
                        float v3 = im(x+3, y, t, c);
                        im(x+2, y, t, c) = DAUB2 * v0 + DAUB1 * v1 + DAUB0 * v2 + DAUB3 * v3;
                        im(x+3, y, t, c) = DAUB3 * v0 - DAUB0 * v1 + DAUB1 * v2 - DAUB2 * v3;
                    }


                    // special case the first two elements using rotation
                    float v0 = saved2nd;
                    float v1 = saved1st;
                    float v2 = im(0, y, t, c);
                    float v3 = im(1, y, t, c);
                    im(0, y, t, c) = DAUB2 * v0 + DAUB1 * v1 + DAUB0 * v2 + DAUB3 * v3;
                    im(1, y, t, c) = DAUB3 * v0 - DAUB0 * v1 + DAUB1 * v2 - DAUB2 * v3;
                }
            }
        }
    }

    // transform in y
    for (int height = 2; height <= im.height; height *= 2) {
        // Collect averages and differences
        Interleave::apply(im.region(0, 0, 0, 0,
                                    im.width, height, im.frames, im.channels),
                          1, 2, 1);

        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int x = 0; x < im.width; x++) {
                    float saved1st = im(x, height-1, t, c);
                    float saved2nd = im(x, height-2, t, c);

                    for (int y = height-4; y >= 0; y-=2) {
                        float v0 = im(x, y, t, c);
                        float v1 = im(x, y+1, t, c);
                        float v2 = im(x, y+2, t, c);
                        float v3 = im(x, y+3, t, c);
                        im(x, y+2, t, c) = DAUB2 * v0 + DAUB1 * v1 + DAUB0 * v2 + DAUB3 * v3;
                        im(x, y+3, t, c) = DAUB3 * v0 - DAUB0 * v1 + DAUB1 * v2 - DAUB2 * v3;
                    }

                    // special case the first two elements using rotation
                    float v0 = saved2nd;
                    float v1 = saved1st;
                    float v2 = im(x, 0, t, c);
                    float v3 = im(x, 1, t, c);
                    im(x, 0, t, c) = DAUB2 * v0 + DAUB1 * v1 + DAUB0 * v2 + DAUB3 * v3;
                    im(x, 1, t, c) = DAUB3 * v0 - DAUB0 * v1 + DAUB1 * v2 - DAUB2 * v3;
                }
            }
        }
    }
}



#include "footer.h"
