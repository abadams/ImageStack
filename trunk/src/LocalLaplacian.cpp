#include "main.h"
#include "LocalLaplacian.h"
#include "Geometry.h"
#include "Arithmetic.h"
#include "Statistics.h"
#include "Filter.h"
#include "File.h"
#include "Func.h"
#include "header.h"

/*
Func LocalLaplacian::pyramidDown(Func im) {
    Func dx = subsampleX(im, 2, -1) + 3*subsampleX(im, 2, 0) + 3*subsampleX(im, 2, 1) + subsampleX(im, 2, 2);
    Func dy = subsampleY(dx, 2, -1) + 3*subsampleY(dx, 2, 0) + 3*subsampleY(dx, 2, 1) + subsampleY(dx, 2, 2);
    return dy;
}

Func LocalLaplacian::pyramidUp(Func im) {
    Func upX = interleaveX(im, im);
    upX = shiftX(upX, -1) + 2*upX + shiftX(upX, 1);
    Func upY = interleaveY(upX, upX);
    upX = shiftY(upY, -1) + 2*upY + shiftY(upY, 1);    
    return upY;
}
*/

Image pyramidDown(Image im) {
    auto zb = Expr::zeroBoundary(im);
    Expr::Func dx = (subsampleX(zb, 2, -1) + 
                     3*subsampleX(zb, 2, 0) + 
                     3*subsampleX(zb, 2, 1) + 
                     subsampleX(zb, 2, 2));
    auto dy = (subsampleY(dx, 2, -1) + 
               3*subsampleY(dx, 2, 0) + 
               3*subsampleY(dx, 2, 1) + 
               subsampleY(dx, 2, 2));
    Image small(im.width/2, im.height/2, im.frames, im.channels);
    small.set(dy/64);
    return small;
    /*
    Image blurry = im.copy();
    FastBlur::apply(blurry, 2, 2, 2);
    return Subsample::apply(blurry, 2, 2, 0, 0);
    */
}

Expr::Func pyramidUp(Image im) {
    //Image larger(w, h, f, im.channels);

    /*
    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < f; t++) {
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    larger(x, y, t, c) = im(x/2, y/2, t/2, c);
                }
            }
        }
    }

    FastBlur::apply(larger, 2, 2, 2);
    */

    auto zb = Expr::zeroBoundary(im);
    Expr::Func upx = interleaveX(shiftX(zb, 1) + 3*zb, 3*zb + shiftX(zb, -1));
    return interleaveY(shiftY(upx, 1) + 3*upx, 3*upx + shiftY(upx, -1))/16;
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

    // Compute a discretized set of K intensities that span the values in the image
    Stats s = Stats(im);
    float target[K];
    for (int i = 0; i < K; i++) {
        target[i] = ((float)i)/(K-1) * (s.maximum() - s.minimum()) + s.minimum();
    }

    float sigma = 1.0f/(target[1] - target[0]);
    alpha /= K;

    // Convert to grayscale
    Image gray(im.width, im.height, im.frames, 1);
    for (int c = 0; c < im.channels; c++) {
        gray += im.channel(c) / im.channels;
    }

    // Compute a Gaussian and Laplacian pyramid for the input
    printf("Computing Gaussian pyramid for input\n");
    Image imPyramid[J];
    Image imLPyramid[J];
    imPyramid[0] = gray;
    for (int j = 1; j < J; j++) {
        imPyramid[j] = pyramidDown(imPyramid[j-1]);
        imLPyramid[j-1] = imPyramid[j-1] - pyramidUp(imPyramid[j]);
    }
    imLPyramid[J-1] = imPyramid[J-1];

    // Make a lookup table for remapping
    // It's the derivative of a Gaussian centered at 1024 with std.dev 256
    Image remap(16*256, 1, 1, 1);
    auto fx = (Expr::X()-8*256) / 256.0f;
    remap.set(alpha*fx*exp(-fx*fx/2.0f));
    Save::apply(remap, "remap.tmp");
    
    // Make a set of K processed images
    printf("Computing different processed images\n");
    Image pyramid[K][J];
    for (int i = 0; i < K; i++) {
        Image p(gray.width, gray.height, gray.frames, 1);

        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width; x++) {
                    float luminance = gray(x, y, t, 0);
                    float v = (luminance - target[i]) * sigma;
                    int idx = clamp((int)(v * 256 + 8*256), 0, remap.width-1);
                    float adjustment = remap(idx, 0); //alpha * v * fastexp(sigma*v*v);
                    p(x, y, t, 0) = luminance + adjustment;
                }
            }
        }

        // Compute a J-level laplacian pyramid per processed image
        pyramid[i][0] = p;
        for (int j = 1; j < J; j++) {
            pyramid[i][j] = pyramidDown(pyramid[i][j-1]);
            pyramid[i][j-1] -= pyramidUp(pyramid[i][j]);
        }
    }

    // Now construct output laplacian pyramid by looking up the
    // Laplacian pyramids as a function of intensity found in the
    // Gaussian pyramid
    printf("Computing laplacian pyramid for output\n");
    Image output(imPyramid[J-1].width, imPyramid[J-1].height, imPyramid[J-1].frames, 1);
    for (int j = J-1; j >= 0; j--) {

        float scale;
        if (beta < 0) {
            // when j = 0, scale =
            scale = ((float)j/(J-1))*(-beta) + 1-(-beta);
        } else {
            scale = (1.0f - ((float)j/(J-1)))*beta + 1-beta;
        }

        printf("%d %f\n", j, scale);

        // Upsample what we have so far
        if (j != J-1) {            
            Image newOutput = Image(imPyramid[j].width, imPyramid[j].height, 
                                    imPyramid[j].frames, imPyramid[j].channels);
            newOutput.set(pyramidUp(output));
            output = newOutput;
        }

        printf("output size: %d %d %d %d\n", output.width, output.height, output.frames, output.channels);

        // Add in the next pyramid level
        for (int t = 0; t < imPyramid[j].frames; t++) {
            for (int y = 0; y < imPyramid[j].height; y++) {
                for (int x = 0; x < imPyramid[j].width; x++) {
                    float luminance = imPyramid[j](x, y, t, 0);
                    
                    luminance -= s.minimum();
                    luminance /= s.maximum() - s.minimum();
                    luminance *= K-1;
                    
                    int K0 = int(luminance);
                    if (K0 < 0) K0 = 0;
                    if (K0 >= K-1) K0 = K-1;
                    
                    int K1 = K0+1;
                    
                    float interp = luminance - K0;
                    
                    float modified =
                        interp * pyramid[K1][j](x, y, t, 0) +
                        (1 - interp) * pyramid[K0][j](x, y, t, 0);
                    
                    output(x, y, t, 0) += 
                      (1-scale) * imLPyramid[j](x, y, t, 0) + scale * modified;
                }
            }
        }

    }

    Save::apply(output, "output.tmp"); 

    // Reintroduce color
    for (int c = 0; c < im.channels; c++) {
        im.channel(c).set((im.channel(c) / gray) * output);
    }
}

#include "footer.h"
