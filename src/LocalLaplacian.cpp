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

Func pyramidDown(Func im) {
    Func dx = (subsampleX(im, 2, -1) + 3*subsampleX(im, 2, 0) + 3*subsampleX(im, 2, 1) + subsampleX(im, 2, 2))/8;
    Func dy = (subsampleY(dx, 2, -1) + 3*subsampleY(dx, 2, 0) + 3*subsampleY(dx, 2, 1) + subsampleY(dx, 2, 2))/8;
    dx.setName("dx");
    dy.setName("dy");
    return dy;
}

Func pyramidUp(Func im) {
    Func upx = interleaveX(3*im + shiftX(im, 1), 3*im + shiftX(im, -1));
    Func upy = interleaveY(3*upx + shiftY(upx, 1), 3*upx + shiftY(upx, -1))/16;
    upx.setName("upx");
    upy.setName("upy");
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
    Func gray = zeroBoundary(im.channel(0) + im.channel(1) + im.channel(2))/3;
    gray.setName("gray");

    // Compute a Gaussian and Laplacian pyramid for the input
    printf("Computing Gaussian pyramid for input\n");
    Func imPyramid[J];
    Func imLPyramid[J];
    imPyramid[0] = gray;
    for (int j = 1; j < J; j++) {
        imPyramid[j] = pyramidDown(imPyramid[j-1]);
        imLPyramid[j-1] = imPyramid[j-1] - pyramidUp(imPyramid[j]);
        
        {
            std::stringstream ss;
            ss << "imPyramid_" << j;
            imPyramid[j].setName(ss.str());
        }
        {
            std::stringstream ss;
            ss << "imLPyramid_" << (j-1);
            imLPyramid[j-1].setName(ss.str());
        }
        imLPyramid[J-1] = imPyramid[J-1];
    }


    // Make a lookup table for remapping
    // It's the derivative of a Gaussian centered at 1024 with std.dev 256
    Image remap(16*256, 1, 1, 1);
    auto fx = (Expr::X()-8*256) / 256.0f;
    remap.set(alpha*fx*exp(-fx*fx/2.0f));
    
    // Make a set of K processed images
    printf("Computing different processed images\n");
    Func pyramid[J];
    Expr::X x; Expr::Y y; Expr::C c;    
    auto diff = (gray(x, y) - minIntensity) / intensityDelta;
    auto idx = clamp(toInt(diff * 256) - c*256 + remap.width/2, 0, remap.width-1);
    pyramid[0] = gray(x, y) + remap(idx);

    // Compute a J-level Gaussian pyramid of the processed images
    for (int j = 1; j < J; j++) {
        pyramid[j] = pyramidDown(pyramid[j-1]);
        std::stringstream ss;
        ss << "pyramid_" << j;
        pyramid[j].setName(ss.str());
    }

    // Compute a J-level laplacian pyramid of the processed images
    Func lPyramid[J];
    lPyramid[J-1] = pyramid[J-1];
    for (int j = J-2; j >= 0; j--) {
        lPyramid[j] = pyramid[j] - pyramidUp(pyramid[j+1]);
        std::stringstream ss;
        ss << "lPyramid_" << j;
        lPyramid[j].setName(ss.str());
        lPyramid[j].eager(); 
    }

    // Now construct output laplacian pyramid by looking up the
    // Laplacian pyramids as a function of intensity found in the
    // Gaussian pyramid
    printf("Computing laplacian pyramid for output\n");
    Func output[J];
    for (int j = J-1; j >= 0; j--) {

        float scale;
        if (beta < 0) {
            // when j = 0, scale =
            scale = ((float)j/(J-1))*(-beta) + 1-(-beta);
        } else {
            scale = (1.0f - ((float)j/(J-1)))*beta + 1-beta;
        }

        printf("%d %f\n", j, scale);

        auto level = (imPyramid[j] - minIntensity)/intensityDelta;
        auto intLevel = clamp(toInt(level), 0, K-2);
        auto interp = level - toFloat(intLevel);
        // Add in the next pyramid level        
        auto modified = 
            interp * lPyramid[j](x, y, intLevel+1) + 
            (1 - interp) * lPyramid[j](x, y, intLevel);
        auto newVal = (1-scale) * imLPyramid[j] + scale * modified;
        if (j == J-1) {
            output[j] = newVal;
        } else {
            output[j] = pyramidUp(output[j+1]) + newVal;
        }

        std::stringstream ss;
        ss << "output_" << j;
        output[j].setName(ss.str());

    }

    // Reintroduce color
    im *= output[0](x, y) / gray(x, y);
}

#include "footer.h"
