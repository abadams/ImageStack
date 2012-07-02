#include "main.h"
#include "LightField.h"
#include "Geometry.h"
#include "Stack.h"
#include "Convolve.h"
#include "Arithmetic.h"
#include "Color.h"
#include "File.h"
#include "Wavelet.h"
#include "Filter.h"
#include "LinearAlgebra.h"
#include "header.h"

void LFFocalStack::help() {
    printf("\n-lffocalstack turns a 4d light field into a 3d focal stack. The five arguments\n"
           "are the lenslet width, height, the minimum alpha, the maximum alpha, and the\n"
           "step size between adjacent depths (alpha is slope in line space).\n\n"
           "Usage: ImageStack -load lf.exr -lffocalstack 16 16 -1 1 0.1 -display\n\n");
}

bool LFFocalStack::test() {
    Image im(256, 256, 1, 1);
    // 16 x 16 x 16 x 16 lightfield
    LightField point(im, 16, 16);
    LFPoint::apply(point, 7, 7, 0.5);
    Image stack = LFFocalStack::apply(point, -1, 1, 0.5);

    // Should be a focused spot at one particular frame
    float spot = stack(7, 7, 3, 0);

    // Zero elsewhere in that frame
    float zero = stack(10, 10, 3, 0);

    // At another frame it should be blurred out
    float gray = stack(10, 10, 0, 0);

    return spot > 0.75 && nearlyEqual(zero, 0) && gray < spot/16 && gray > spot/256;
}

void LFFocalStack::parse(vector<string> args) {
    assert(args.size() == 5, "-lffocalstack takes five arguments.\n");
    LightField lf(stack(0), readInt(args[0]), readInt(args[1]));
    Image im = apply(lf, readFloat(args[2]), readFloat(args[3]), readFloat(args[4]));
    pop();
    push(im);
}


Image LFFocalStack::apply(LightField lf, float minAlpha, float maxAlpha, float deltaAlpha) {
    assert(lf.image.frames == 1, "Can only turn a single light field into a focal stack\n");

    int frames = 0;
    for (float alpha = minAlpha; alpha <= maxAlpha; alpha += deltaAlpha) { frames++; }

    Image out(lf.xSize, lf.ySize, frames, lf.image.channels);

    Image view(lf.xSize, lf.ySize, 1, lf.image.channels);

    int t = 0;
    for (float alpha = minAlpha; alpha <= maxAlpha; alpha += deltaAlpha) {
        printf("computing frame %i\n", t+1);

        // Extract, shift, and accumulate each view
        for (int v = 0; v < lf.vSize; v++) {
            for (int u = 0; u < lf.uSize; u++) {
                // Get the view
                for (int y = 0; y < lf.ySize; y++) {
                    for (int x = 0; x < lf.xSize; x++) {
                        for (int c = 0; c < lf.image.channels; c++) {
                            view(x, y, c) = lf(x, y, u, v, c);
                        }
                    }
                }

                // Shift it
                if (alpha != 0) {
                    view = Translate::apply(view, -(u-(lf.uSize-1)*0.5)*alpha, -(v-(lf.vSize-1)*0.5)*alpha);
                }

                // Accumulate it into the output
                out.frame(t) += view;
            }
        }

        // Filter if necessary
        if (fabs(alpha) > 1) {
            out.frame(t).set(LanczosBlur::apply(out.frame(t), fabs(alpha), fabs(alpha), 0));
        }
        t++;
    }

    // renormalize
    out /= lf.uSize * lf.vSize;

    return out;
}

void LFPoint::help() {
    printf("\n-lfpoint colors a single 3d point white in the given light field. The five\n"
           "arguments are the light field u, v, resolution, and then the x, y, and z\n"
           "coordinates of the point. x and y should be in the range [0, 1], while z\n"
           "is disparity. z = 0 will be at the focal plane.\n\n"
           "Usage: ImageStack -load lf.exr -lfpoint 16 16 0.5 0.5 0.1 -save newlf.exr\n\n");
}

bool LFPoint::test() {
    // tested by LFFocalStack
    return true;
}

void LFPoint::parse(vector<string> args) {
    assert(args.size() == 5, "-lfpoint takes five arguments\n");
    LightField lf(stack(0), readInt(args[0]), readInt(args[1]));
    apply(lf, readFloat(args[2]), readFloat(args[3]), readFloat(args[4]));
}

void LFPoint::apply(LightField lf, float px, float py, float pz) {
    for (int v = 0; v < lf.vSize; v++) {
        for (int u = 0; u < lf.uSize; u++) {
            float pu = u - (lf.uSize-1) * 0.5;
            float pv = v - (lf.vSize-1) * 0.5;
            // figure out the correct x y
            int x = (int)floorf(px + pz*pu + 0.5);
            int y = (int)floorf(py + pz*pv + 0.5);
            if (x < 0 || x >= lf.xSize || y < 0 || y >= lf.ySize) continue;
            for (int c = 0; c < lf.image.channels; c++) {
                lf(x, y, u, v, c) = 1;
            }
        }
    }
}

#include "footer.h"
