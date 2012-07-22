#include "main.h"
#include "LocalLaplacian.h"
#include "Geometry.h"
#include "Arithmetic.h"
#include "Statistics.h"
#include "Filter.h"
#include "File.h"
#include "Func.h"
#include "header.h"

using namespace Expr;

Image pyramidDown(Image im) {
    Image dx = (subsampleX(zeroBoundary(im), 2, -1) + 3*subsampleX(im, 2, 0) + 
                3*subsampleX(zeroBoundary(im), 2, 1) + subsampleX(zeroBoundary(im), 2, 2));
    Image dy = (subsampleY(zeroBoundary(dx), 2, -1) + 3*subsampleY(dx, 2, 0) + 
                3*subsampleY(zeroBoundary(dx), 2, 1) + subsampleY(zeroBoundary(dx), 2, 2))/64;
    return dy;
}

Image pyramidUp(Image im) {
    Image upx(im.width*2+1, im.height, im.frames, im.channels);
    Image upy(upx.width, im.height*2+1, im.frames, im.channels);
    upx.set(interleaveX(3*zeroBoundary(im) + shiftX(zeroBoundary(im), 1), 
                        3*zeroBoundary(im) + shiftX(zeroBoundary(im), -1)));
    upy.set(interleaveY(3*zeroBoundary(upx) + shiftY(zeroBoundary(upx), 1), 
                        3*zeroBoundary(upx) + shiftY(zeroBoundary(upx), -1))/16);
    return upy;
}

void LocalLaplacian::help() {
    pprintf("-locallaplacian modifies contrast at various scales. It is similar to"
            " the clarity slider in Photoshop. This operation is an implementation"
            " of Fast and Robust Pyramid-based Image Processing by Aubry et al,"
            " which is an acceleration of Local Laplacian Filters by Paris et"
            " al.\n\n"
            "The first argument specifies how much of an effect to apply. 0"
            " gives no effect, 1 produces a moderate amount of contrast"
            " enhancement, and -1 produces a piecewise flattening. The second"
            " argument specifies how the effect should change with respect to"
            " scale. 0 applies the same effect at all scales, -1 applies the effect"
            " only at coarse scales, and 1 applies the effect only at fine"
            " scales. Values larger than one actually reverse the effect at fine or"
            " coarse scales. In fact -locallaplacian 1 2 is an effective"
            " tone-mapper because it amplifies contrast at fine scales and reduces"
            " it at coarse scales.\n"
            "\n"
            "Usage: ImageStack -load input.jpg -locallaplacian 1 0 -save boosted.jpg\n");
}

bool LocalLaplacian::test() {
    Image im = Downsample::apply(Load::apply("pics/dog1.jpg"), 2, 2, 1);
    Stats si(im);
    si.variance();
    LocalLaplacian::apply(im, 1.2, 0.2);
    Stats se(im);
    return (se.minimum() < si.minimum() &&
            se.maximum() > si.maximum() &&
            se.variance() > si.variance());
}

void LocalLaplacian::parse(vector<string> args) {
    assert(args.size() == 2, "-locallaplacian takes two arguments");
    Image im = stack(0);
    apply(stack(0), readFloat(args[0]), readFloat(args[1]));
}

void LocalLaplacian::apply(Image im, float alpha, float beta) {
    const int K = 8, J = 8;

    assert(im.channels == 3, "-locallaplacian only works on three-channel images\n");

    // For multi-frame images, process each frame independently
    if (im.frames > 1) {
        for (int t = 0; t < im.frames; t++) {
            apply(im.frame(t), alpha, beta);
        }
        return;
    }

    // Compute a discretized set of K intensities that span the values in the image
    Stats s = Stats(im);
    float minIntensity = s.minimum();
    float intensityDelta = (s.maximum() - s.minimum()) / (K-1);
    alpha /= (K-1);

    // Convert to grayscale and add a boundary condition
    Image gray = (im.channel(0) + im.channel(1) + im.channel(2))/3;

    // Compute a Gaussian and Laplacian pyramid for the input
    Image imPyramid[J];
    Image imLPyramid[J];
    imPyramid[0] = gray;
    for (int j = 1; j < J; j++) {
        imPyramid[j] = pyramidDown(imPyramid[j-1]);
        imLPyramid[j-1] = imPyramid[j-1] - zeroBoundary(pyramidUp(imPyramid[j]));
    }
    imLPyramid[J-1] = imPyramid[J-1];

    // Make a lookup table for remapping
    // It's the derivative of a Gaussian centered at 1024 with std.dev 256
    Image remap(16*256, 1, 1, 1);
    auto fx = (Expr::X()-8*256) / 256.0f;
    remap.set(alpha*fx*exp(-fx*fx/2.0f));
    
    // Make a set of K processed images
    Image pyramid[J];
    pyramid[0] = Image(gray.width, gray.height, 1, K);
    Expr::X x; Expr::Y y; Expr::C c;
    auto diff = (gray(x, y) - minIntensity) / intensityDelta;
    auto idx = clamp(toInt(diff * 256) - c*256 + remap.width/2, 0, remap.width-1);        
    pyramid[0].set(gray(x, y) + remap(idx));

    // Compute a J-level Laplacian pyramid of the processed images
    for (int j = 1; j < J; j++) {
        pyramid[j] = pyramidDown(pyramid[j-1]);
        pyramid[j-1] = pyramid[j-1] - zeroBoundary(pyramidUp(pyramid[j]));
    }

    // Now construct output laplacian pyramid by looking up the
    // Laplacian pyramids as a function of intensity found in the
    // Gaussian pyramid
    Image output;
    for (int j = J-1; j >= 0; j--) {

        float scale;
        if (beta < 0) {
            scale = ((float)j/(J-1))*(-beta) + 1-(-beta);
        } else {
            scale = (1.0f - ((float)j/(J-1)))*beta + 1-beta;
        }

        auto level = (imPyramid[j] - minIntensity)/intensityDelta;
        auto intLevel = clamp(toInt(level), 0, K-2);
        auto interp = level - toFloat(intLevel);
        // Add in the next pyramid level        
        auto modified = 
            interp * pyramid[j](x, y, intLevel+1) + 
            (1 - interp) * pyramid[j](x, y, intLevel);
        auto newVal = (1-scale) * imLPyramid[j] + scale * modified;
        if (j == J-1) {
            output = newVal;
        } else {
            output = zeroBoundary(pyramidUp(output)) + newVal;
        }
    }

    // Reintroduce color
    output /= gray;
    im *= output(x, y, 0);
}

#include "footer.h"
