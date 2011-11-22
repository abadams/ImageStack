#include "main.h"
#include "LocalLaplacian.h"
#include "Geometry.h"
#include "Arithmetic.h"
#include "Statistics.h"
#include "Filter.h"
#include "File.h"
#include "header.h"

Image LocalLaplacian::pyramidDown(Window im) {
    Image smaller((im.width+1)/2, (im.height+1)/2, (im.frames+1)/2, im.channels);

    Image blurry(im);
    FastBlur::apply(blurry, 1, 1, 1);

    for (int t = 0; t < im.frames; t++) {
        for (int y = 0; y < im.height; y++) {
            for (int x = 0; x < im.width; x++) {
                for (int c = 0; c < im.channels; c++) {
                    smaller(x/2, y/2, t/2)[c] += blurry(x, y, t)[c];
                }
            }
        }
    }

    float scale[4] = {0.125f, 0.25f, 0.5f, 1.0f};

    for (int t = 0; t < smaller.frames; t++) {
        int tFactor = ((t*2+1) < im.frames) ? 0 : 1;
        for (int y = 0; y < smaller.height; y++) {
            int yFactor = ((y*2+1) < im.height) ? tFactor : (tFactor+1);
            for (int x = 0; x < smaller.width; x++) {
                int xFactor = ((x*2+1) < im.width) ? (yFactor) : (yFactor+1);
                for (int c = 0; c < smaller.channels; c++) {
                    smaller(x, y, t)[c] *= scale[xFactor];
                }
            }
        }
    }
    
    return smaller;
}

Image LocalLaplacian::pyramidUp(Window im, int w, int h, int f) {
    Image larger(w, h, f, im.channels);

    for (int t = 0; t < f; t++) {
        for (int y = 0; y < h; y++) {
            for (int x = 0; x < w; x++) {
                for (int c = 0; c < im.channels; c++) {
                    larger(x, y, t)[c] = im(x/2, y/2, t/2)[c];
                }
            }
        }
    }

    FastBlur::apply(larger, 1, 1, 1, false);

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

void LocalLaplacian::parse(vector<string> args) {
    assert(args.size() == 2, "-locallaplacian takes two arguments");
    Image im = stack(0);
    pop();
    push(apply(im, readFloat(args[0]), readFloat(args[1]))); 
}

Image LocalLaplacian::apply(Window im, float alpha, float beta) {
    const int K = 8, J = 8;

    // Compute a discretized set of K intensities that span the values in the image
    Stats s = Stats(im);
    float target[K];
    for (int i = 0; i < K; i++) {
        target[i] = ((float)i)/(K-1) * (s.maximum() - s.minimum()) + s.minimum();
    }

    float sigma = 1.0f/(target[1] - target[0]);
    sigma = - sigma * sigma * 0.5f;

    // Make a set of K processed images
    printf("Computing different processed images\n");
    Image processed[K];
    for (int i = 0; i < K; i++) {
        Image p(im);
        
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width; x++) {
                    float luminance = 0;
                    for (int c = 0; c < im.channels; c++) {
                        luminance += p(x, y, t)[c];
                    }
                    luminance /= im.channels;

                    // Attract (or repel for negative alpha) the luminance towards the target value
                    //float newLuminance = alpha * luminance + (1 - alpha) * target[i];

                    float v = luminance - target[i];
                    float adjustment = alpha * v * fastexp(sigma*v*v);

                    //printf("At level %d, %f -> %f\n", i, luminance, luminance + adjustment);

                    for (int c = 0; c < im.channels; c++) {
                        p(x, y, t)[c] += adjustment;
                    }
                }
            }
        }
        processed[i] = p;
    }

    // Compute a J-level laplacian pyramid per processed image
    printf("Computing their laplacian pyramids\n");
    Image pyramid[K][J];
    for (int i = 0; i < K; i++) {
        pyramid[i][0] = processed[i];
        for (int j = 1; j < J; j++) {
            int oldW = pyramid[i][j-1].width;
            int oldH = pyramid[i][j-1].height;
            int oldF = pyramid[i][j-1].frames;
            pyramid[i][j] = pyramidDown(pyramid[i][j-1]);
            Subtract::apply(pyramid[i][j-1], pyramidUp(pyramid[i][j], oldW, oldH, oldF));
        }
    }

    // Now compute a Gaussian and Laplacian pyramid for the input
    printf("Computing Gaussian pyramid for input\n");
    Image imPyramid[J];
    Image imLPyramid[J];
    imPyramid[0] = Image(im);
    imLPyramid[0] = Image(im);
    for (int j = 1; j < J; j++) {
        int oldW = imPyramid[j-1].width;
        int oldH = imPyramid[j-1].height;
        int oldF = imPyramid[j-1].frames;
        imPyramid[j] = pyramidDown(imPyramid[j-1]);
        imLPyramid[j] = imPyramid[j].copy();
        Subtract::apply(imLPyramid[j-1], pyramidUp(imPyramid[j], oldW, oldH, oldF));
    }

    // Now construct output laplacian pyramid by looking up the
    // Laplacian pyramids as a function of intensity found in the
    // Gaussian pyramid
    printf("Computing laplacian pyramid for input\n");
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
                        luminance += imPyramid[j](x, y, t)[c];
                    }
                    luminance /= im.channels;

                    luminance -= s.minimum();
                    luminance /= s.maximum() - s.minimum();
                    luminance *= K-1;
                    
                    int K0 = int(luminance);
                    if (K0 < 0) K0 = 0;
                    if (K0 >= K-1) K0 = K-1;

                    int K1 = K0+1;

                    float alpha = luminance - K0;

                    for (int c = 0; c < im.channels; c++) {
                        float modified = 
                            alpha * pyramid[K1][j](x, y, t)[c] +
                            (1 - alpha) * pyramid[K0][j](x, y, t)[c];

                        imLPyramid[j](x, y, t)[c] = 
                            (1-scale) * imLPyramid[j](x, y, t)[c] + scale * modified;

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
        Add::apply(output, imLPyramid[j]);
    }

    return output;
}

#include "footer.h"
