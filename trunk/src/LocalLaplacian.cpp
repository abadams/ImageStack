#include "main.h"
#include "LocalLaplacian.h"
#include "Geometry.h"
#include "Arithmetic.h"
#include "Statistics.h"
#include "Filter.h"
#include "File.h"
#include "header.h"

/*
Func LocalLaplacian::pyramidDown(Func im) {
    Func dx = subsampleX(im, 2, -1) + 3*subsampleX(im, 2, 0) + 3*subsampleX(im, 2, 1) + subsampleX(im, 2, 2);
    Func dy = subsampleY(dy, 2, -1) + 3*subsampleY(dy, 2, 0) + 3*subsampleY(dy, 2, 1) + subsampleY(dy, 2, 2);
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

Image LocalLaplacian::pyramidDown(Image im) {
    Image blurry = im.copy();
    FastBlur::apply(blurry, 2, 2, 2);
    return Subsample::apply(blurry, 2, 2, 0, 0);
}

Image LocalLaplacian::pyramidUp(Image im, int w, int h, int f) {
    Image larger(w, h, f, im.channels);

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

    return larger;
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
    Image enhance = LocalLaplacian::apply(im, 1.2, 0.2);
    Stats si(im), se(enhance);
    return (se.minimum() < si.minimum() &&
            se.maximum() > si.maximum() &&
            se.variance() > si.variance() &&
            nearlyEqual(si.mean(), se.mean()));
}

void LocalLaplacian::parse(vector<string> args) {
    assert(args.size() == 2, "-locallaplacian takes two arguments");
    Image im = stack(0);
    pop();
    push(apply(im, readFloat(args[0]), readFloat(args[1])));
}

Image LocalLaplacian::apply(Image im, float alpha, float beta) {
    const int K = 8, J = 8;

    // Compute a discretized set of K intensities that span the values in the image
    Stats s = Stats(im);
    float target[K];
    for (int i = 0; i < K; i++) {
        target[i] = ((float)i)/(K-1) * (s.maximum() - s.minimum()) + s.minimum();
    }

    float sigma = 1.0f/(target[1] - target[0]);
    sigma = - sigma * sigma * 0.5f;

    // Compute a Gaussian and Laplacian pyramid for the input
    printf("Computing Gaussian pyramid for input\n");
    Image imPyramid[J];
    Image imLPyramid[J];
    imPyramid[0] = im;
    imLPyramid[0] = im.copy();
    for (int j = 1; j < J; j++) {
        int oldW = imPyramid[j-1].width;
        int oldH = imPyramid[j-1].height;
        int oldF = imPyramid[j-1].frames;
        imPyramid[j] = pyramidDown(imPyramid[j-1]);
        imLPyramid[j] = imPyramid[j].copy();
        imLPyramid[j-1] -= pyramidUp(imPyramid[j], oldW, oldH, oldF);
    }

    // Make a set of K processed images
    printf("Computing different processed images\n");
    Image processed[K];
    Image pyramid[K][J];
    for (int i = 0; i < K; i++) {
        Image p = im.copy();


        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width; x++) {
                    float luminance = 0;
                    for (int c = 0; c < im.channels; c++) {
                        luminance += p(x, y, t, c);
                    }
                    luminance /= im.channels;

                    // Attract (or repel for negative alpha) the luminance towards the target value
                    // float newLuminance = alpha * luminance + (1 - alpha) * target[i];

                    float v = luminance - target[i];
                    float adjustment = alpha * v * fastexp(sigma*v*v);

                    //printf("At level %d, %f -> %f\n", i, luminance, luminance + adjustment);

                    for (int c = 0; c < im.channels; c++) {
                        p(x, y, t, c) += adjustment;
                    }
                }
            }
        }
        processed[i] = p;

        // Compute a J-level laplacian pyramid per processed image
        pyramid[i][0] = processed[i];
        for (int j = 1; j < J; j++) {
            int oldW = pyramid[i][j-1].width;
            int oldH = pyramid[i][j-1].height;
            int oldF = pyramid[i][j-1].frames;
            pyramid[i][j] = pyramidDown(pyramid[i][j-1]);
            pyramid[i][j-1] -= pyramidUp(pyramid[i][j], oldW, oldH, oldF);
        }
    }

    // Now construct output laplacian pyramid by looking up the
    // Laplacian pyramids as a function of intensity found in the
    // Gaussian pyramid
    printf("Computing laplacian pyramid for output\n");
    for (int j = 0; j < J; j++) {

        float scale;
        if (beta < 0) {
            // when j = 0, scale =
            scale = ((float)j/(J-1))*(-beta) + 1-(-beta);
        } else {
            scale = (1.0f - ((float)j/(J-1)))*beta + 1-beta;
        }

        printf("%d %f\n", j, scale);

        for (int t = 0; t < imPyramid[j].frames; t++) {
            for (int y = 0; y < imPyramid[j].height; y++) {
                for (int x = 0; x < imPyramid[j].width; x++) {
                    float luminance = 0;
                    for (int c = 0; c < im.channels; c++) {
                        luminance += imPyramid[j](x, y, t, c);
                    }
                    luminance /= im.channels;

                    luminance -= s.minimum();
                    luminance /= s.maximum() - s.minimum();
                    luminance *= K-1;

                    int K0 = int(luminance);
                    if (K0 < 0) K0 = 0;
                    if (K0 >= K-1) K0 = K-1;

                    int K1 = K0+1;

                    float interp = luminance - K0;

                    for (int c = 0; c < im.channels; c++) {
                        float modified =
                            interp * pyramid[K1][j](x, y, t, c) +
                            (1 - interp) * pyramid[K0][j](x, y, t, c);

                        imLPyramid[j](x, y, t, c) =
                            (1-scale) * imLPyramid[j](x, y, t, c) + scale * modified;

                    }
                }
            }
        }
    }

    // Now collapse the output laplacian pyramid
    printf("Collapsing laplacian pyramid down to output image\n");
    Image output = imLPyramid[J-1];
    for (int j = J-2; j >= 0; j--) {
        output = pyramidUp(output, imLPyramid[j].width, imLPyramid[j].height, imLPyramid[j].frames);
        output += imLPyramid[j];
    }

    return output;
}

#include "footer.h"
