
// CS448F final project
// Implementation of PatchMatch algorithm and its applications
// Sung Hee Park (shpark7@stanford.edu)


#include "main.h"
#include "File.h"
#include "Geometry.h"
#include "PatchMatch.h"
#include "Arithmetic.h"
#include "Calculus.h"
#include "Statistics.h"
#include "Filter.h"
#include "Paint.h"
#include "Prediction.h"
#include "Display.h"
#include "header.h"
// PATCHMATCH =============================================================//

void PatchMatch::help() {

    printf("-patchmatch computes approximate nearest neighbor field from the top\n"
           "image on the stack to the second image on the stack, using the\n"
           "algorithm from the PatchMatch SIGGRAPH 2009 paper. This operation\n"
           "requires two input images which may have multiple frames.\n"
           "It returns an image with four channels. First three channels \n"
           "correspond to x, y, t coordinate of closest patch and \n"
           "fourth channels contains the sum of squared differences \n"
           "between patches. \n"
           "\n"
           " arguments [numIter] [patchSize]\n"
           "  - numIter : number of iterations performed. (default: 5)\n"
           "  - patchSize : size of patch. (default: 7, 7x7 square patch)\n"
           " You can omit some arguments from right to use default values.\n"
           "\n"
           "Usage: ImageStack -load target.jpg -load source.jpg -patchmatch -save match.tmp\n\n");
}

bool PatchMatch::test() {
    // Try a trivial example
    Image dog = Downsample::apply(Load::apply("pics/dog1.jpg"), 2, 2, 1);
    Image shifted = Translate::apply(dog, 10, 5, 0);
    Image flow = PatchMatch::apply(dog, shifted, 5, 7);
    Image correct(dog.width-20, dog.height-20, 1, 4);
    correct.setChannels(Expr::X() + 10, Expr::Y() + 5, 0.0f, 0.0f);
    return nearlyEqual(flow.region(0, 0, 0, 0, dog.width-20, dog.height-20, 1, 4), correct);
}

void PatchMatch::parse(vector<string> args) {

    int numIter = 5, patchSize = 7;

    assert(args.size() <= 2, "-patchmatch takes two or fewer arguments\n");
    if (args.size() == 2) {
        numIter = readInt(args[0]);
        patchSize = (int) readInt(args[1]);
    } else if (args.size() == 1) {
        numIter = readInt(args[0]);
    }

    Image result;

    result = apply(stack(0), stack(1), numIter, patchSize);

    push(result);
}

Image PatchMatch::apply(Image source, Image target, int iterations, int patchSize) {
    return apply(source, target, Image(), iterations, patchSize);
}

Image PatchMatch::apply(Image source, Image target, Image mask, int iterations, int patchSize) {

    if (mask.defined()) {
        assert(target.width == mask.width &&
               target.height == mask.height &&
               target.frames == mask.frames,
               "Mask must have the same dimensions as the target\n");
        assert(mask.channels == 1,
               "Mask must have a single channel\n");
    }
    assert(iterations > 0, "Iterations must be a strictly positive integer\n");
    assert(patchSize >= 3 && (patchSize & 1), "Patch size must be at least 3 and odd\n");

    // convert patch diameter to patch radius
    patchSize /= 2;

    // For each source pixel, output a 3-vector to the best match in
    // the target, with an error as the last channel.
    Image out(source.width, source.height, source.frames, 4);

    // Iterate over source frames, finding a match in the target where
    // the mask is high

    for (int t = 0; t < source.frames; t++) {
        // INITIALIZATION - uniform random assignment
        for (int y = 0; y < source.height; y++) {
            for (int x = 0; x < source.width; x++) {
                int dx = randomInt(patchSize, target.width-patchSize-1);
                int dy = randomInt(patchSize, target.height-patchSize-1);
                int dt = randomInt(0, target.frames-1);
                out(x, y, t, 0) = dx;
                out(x, y, t, 1) = dy;
                out(x, y, t, 2) = dt;
                out(x, y, t, 3) = distance(source, target, mask,
                                           x, y, t,
                                           dx, dy, dt,
                                           patchSize, HUGE_VAL);
            }
        }
    }

    bool forwardSearch = true;

    Image dx = out.channel(0), dy = out.channel(1), dt = out.channel(2), error = out.channel(3);

    for (int i = 0; i < iterations; i++) {

        //printf("Iteration %d\n", i);

        // PROPAGATION
        if (forwardSearch) {
            // Forward propagation - compare left, center and up
            for (int t = 0; t < source.frames; t++) {
                for (int y = 1; y < source.height; y++) {
                    for (int x = 1; x < source.width; x++) {
                        if (error(x, y, t, 0) > 0) {
                            float distLeft = distance(source, target, mask,
                                                      x, y, t,
                                                      dx(x-1, y, t, 0)+1,
                                                      dy(x-1, y, t, 0),
                                                      dt(x-1, y, t, 0),
                                                      patchSize, error(x, y, t, 0));

                            if (distLeft < error(x, y, t, 0)) {
                                dx(x, y, t, 0) = dx(x-1, y, t, 0)+1;
                                dy(x, y, t, 0) = dy(x-1, y, t, 0);
                                dt(x, y, t, 0) = dt(x-1, y, t, 0);
                                error(x, y, t, 0) = distLeft;
                            }

                            float distUp = distance(source, target, mask,
                                                    x, y, t,
                                                    dx(x, y-1, t, 0),
                                                    dy(x, y-1, t, 0)+1,
                                                    dt(x, y-1, t, 0),
                                                    patchSize, error(x, y, t, 0));

                            if (distUp < error(x, y, t, 0)) {
                                dx(x, y, t, 0) = dx(x, y-1, t, 0);
                                dy(x, y, t, 0) = dy(x, y-1, t, 0)+1;
                                dt(x, y, t, 0) = dt(x, y-1, t, 0);
                                error(x, y, t, 0) = distUp;
                            }
                        }

                        // TODO: Consider searching across time as well

                    }
                }
            }

        } else {
            // Backward propagation - compare right, center and down
            for (int t = source.frames-1; t >= 0; t--) {
                for (int y = source.height-2; y >= 0; y--) {
                    for (int x = source.width-2; x >= 0; x--) {
                        if (error(x, y, t, 0) > 0) {
                            float distRight = distance(source, target, mask,
                                                       x, y, t,
                                                       dx(x+1, y, t, 0)-1,
                                                       dy(x+1, y, t, 0),
                                                       dt(x+1, y, t, 0),
                                                       patchSize, error(x, y, t, 0));

                            if (distRight < error(x, y, t, 0)) {
                                dx(x, y, t, 0) = dx(x+1, y, t, 0)-1;
                                dy(x, y, t, 0) = dy(x+1, y, t, 0);
                                dt(x, y, t, 0) = dt(x+1, y, t, 0);
                                error(x, y, t, 0) = distRight;
                            }

                            float distDown = distance(source, target, mask,
                                                      x, y, t,
                                                      dx(x, y+1, t, 0),
                                                      dy(x, y+1, t, 0)-1,
                                                      dt(x, y+1, t, 0),
                                                      patchSize, error(x, y, t, 0));

                            if (distDown < error(x, y, t, 0)) {
                                dx(x, y, t, 0) = dx(x, y+1, t, 0);
                                dy(x, y, t, 0) = dy(x, y+1, t, 0)-1;
                                dt(x, y, t, 0) = dt(x, y+1, t, 0);
                                error(x, y, t, 0) = distDown;
                            }
                        }

                        // TODO: Consider searching across time as well

                    }
                }
            }
        }

        forwardSearch = !forwardSearch;

        // RANDOM SEARCH
        for (int t = 0; t < source.frames; t++) {
            for (int y = 0; y < source.height; y++) {
                for (int x = 0; x < source.width; x++) {
                    if (error(x, y, t, 0) > 0) {

                        int radius = target.width > target.height ? target.width : target.height;

                        // search an exponentially smaller window each iteration
                        while (radius > 8) {
                            // Search around current offset vector (distance-weighted)

                            // clamp the search window to the image
                            int minX = (int)dx(x, y, t, 0) - radius;
                            int maxX = (int)dx(x, y, t, 0) + radius + 1;
                            int minY = (int)dy(x, y, t, 0) - radius;
                            int maxY = (int)dy(x, y, t, 0) + radius + 1;
                            if (minX < 0) { minX = 0; }
                            if (maxX > target.width) { maxX = target.width; }
                            if (minY < 0) { minY = 0; }
                            if (maxY > target.height) { maxY = target.height; }

                            int randX = randomInt(minX, maxX-1);
                            int randY = randomInt(minY, maxY-1);
                            int randT = randomInt(0, target.frames - 1);
                            float dist = distance(source, target, mask,
                                                  x, y, t,
                                                  randX, randY, randT,
                                                  patchSize, error(x, y, t, 0));
                            if (dist < error(x, y, t, 0)) {
                                dx(x, y, t, 0) = randX;
                                dy(x, y, t, 0) = randY;
                                dt(x, y, t, 0) = randT;
                                error(x, y, t, 0) = dist;
                            }

                            radius >>= 1;

                        }
                    }
                }
            }
        }
    }

    return out;
}

float PatchMatch::distance(Image source, Image target, Image mask,
                           int sx, int sy, int st,
                           int tx, int ty, int tt,
                           int patchSize, float threshold) {

    // Do not use patches on boundaries
    if (tx < patchSize || tx >= target.width-patchSize ||
        ty < patchSize || ty >= target.height-patchSize) {
        return HUGE_VAL;
    }

    // Compute distance between patches
    // Average L2 distance in RGB space
    float dist = 0;

    int x1 = max(-patchSize, -sx, -tx);
    int x2 = min(patchSize, -sx+source.width-1, -tx+target.width-1);
    int y1 = max(-patchSize, -sy, -ty);
    int y2 = min(patchSize, -sy+source.height-1, -ty+target.height-1);

    for (int c = 0; c < target.channels; c++) {
        for (int y = y1; y <= y2; y++) {
            for (int x = x1; x <= x2; x++) {

                // Don't stray outside the mask
                if (mask.defined() && mask(tx+x, ty+y, tt, 0) < 1) return HUGE_VAL;

                float delta = source(sx+x, sy+y, st, c) - target(tx+x, ty+y, tt, c);
                dist += delta * delta;

                // Early termination
                if (dist > threshold) {return HUGE_VAL;}
            }
        }
    }

    return dist;
}


// BIDIRECTIONAL SIMILARITY =====================================================//


void BidirectionalSimilarity::help() {
    pprintf("-bidirectionalsimilarity reconstructs the top image on the stack using"
            " patches from the second image on the stack, by enforcing coherence"
            " (every patch in the output must look like a patch from the input) and"
            " completeness (every patch from the input must be represented somewhere"
            " in the output). The first argument is a number between zero and one,"
            " which trades off between favoring coherence only (at zero), and"
            " completeness only (at one). It defaults to 0.5. The second arguments"
            " specifies the number of iterations that should be performed, and"
            " defaults to five. Bidirectional similarity uses patchmatch as the"
            " underlying nearest-neighbour-field algorithm, and the third argument"
            " specifies how many iterations of patchmatch should be performed each"
            " time it is run. This also defaults to five.\n"
            "\n"
            "This is an implementation of the paper \"Summarizing visual data using"
            " bidirectional similarity\" by Simakov et al. from CVPR 2008.\n"
            "\n"
            "Usage: ImageStack -load source.jpg -load target.jpg -bidirectional 0.5 -display\n");
}

bool BidirectionalSimilarity::test() {
    // Retarget a duplicated image down to a smaller version of itself
    Image im = Downsample::apply(Load::apply("pics/dog1.jpg"), 8, 8, 1);
    Image source = Adjoin::apply(im, im, 'x');
    source = Adjoin::apply(source, source, 'y');
    Image noisy = im.copy();
    Noise::apply(noisy, -0.3, 0.3);
    // make sure we added enough noise
    if (nearlyEqual(noisy, im)) return false;

    // We should be able to reconstruct a clean dog by resynthesizing it from the source
    BidirectionalSimilarity::apply(source, noisy, Image(), Image(), 0.5, 5, 5);
    return nearlyEqual(noisy, im);
}

void BidirectionalSimilarity::parse(vector<string> args) {

    float alpha = 0.5;
    int numIter = 5;
    int numIterPM = 5;

    assert(args.size() <= 3, "-bidirectional takes three or fewer arguments\n");
    if (args.size() == 3) {
        alpha = readFloat(args[0]);
        numIter = readFloat(args[1]);
        numIterPM = readFloat(args[2]);
    } else if (args.size() == 2) {
        alpha = readFloat(args[0]);
        numIter = readFloat(args[1]);
    } else if (args.size() == 1) {
        alpha = readFloat(args[0]);
    }

    apply(stack(1), stack(0), Image(), Image(), alpha, numIter, numIterPM);
}


// Reconstruct the portion of the target where the mask is high, using
// the portion of the source where its mask is high. Source and target
// masks are allowed to be null Images.
void BidirectionalSimilarity::apply(Image source, Image target,
                                    Image sourceMask, Image targetMask,
                                    float alpha, int numIter, int numIterPM) {


    const int patchSize = 5;

    // TODO: intelligently crop the input to where the mask is high +
    // patch radius on each side

    // Precompute average patch weights
    Image sourceWeight, targetWeight;
    if (sourceMask.defined()) {
        sourceWeight = sourceMask.copy();
        RectFilter::apply(sourceWeight, patchSize, patchSize, 1);
    }
    if (targetMask.defined()) {
        targetWeight = targetMask.copy();
        RectFilter::apply(targetWeight, patchSize, patchSize, 1);
    }

    // recurse
    if (source.width > 32 && source.height > 32 && target.width > 32 && target.height > 32) {
        Image smallSource = Resample::apply(source, source.width/2, source.height/2, source.frames);
        Image smallTarget = Resample::apply(target, target.width/2, target.height/2, target.frames);

        Image smallSourceMask;
        Image smallTargetMask;
        if (sourceMask.defined()) {
            smallSourceMask = Downsample::apply(sourceMask, 2, 2, 1);
        }

        if (targetMask.defined()) {
            smallTargetMask = Downsample::apply(targetMask, 2, 2, 1);
        }

        apply(smallSource, smallTarget, smallSourceMask, smallTargetMask, alpha, numIter, numIterPM);

        Image newTarget = Resample::apply(smallTarget, target.width, target.height, target.frames);

        if (targetMask.defined()) {
            Composite::apply(target, newTarget, targetMask);
        } else {
            target.set(newTarget);
            //newTarget = target.copy();
        }
    }

    printf("%dx%d ", target.width, target.height); fflush(stdout);
    for (int i = 0; i < numIter; i++) {
        printf("."); fflush(stdout);

        // The homogeneous output for this iteration
        Image out(target.width, target.height, target.frames, target.channels+1);

        if (alpha != 0) {

            // COMPLETENESS TERM
            Image completeMatch = PatchMatch::apply(source, target, targetMask, numIterPM, patchSize);

            // For every patch in the source, splat it onto the
            // nearest match in the target, weighted by the source
            // mask and also by the inverse of the patch distance
            for (int t = 0; t < source.frames; t++) {
                for (int y = 0; y < source.height; y++) {
                    for (int x = 0; x < source.width; x++) {

                        float patchWeight = sourceWeight.defined() ? sourceWeight(x, y, t, 0) : 1;

                        // Don't use source patches that aren't completely defined
                        if (patchWeight > 0.99) {

                            int dstX = (int)completeMatch(x, y, t, 0);
                            int dstY = (int)completeMatch(x, y, t, 1);
                            int dstT = (int)completeMatch(x, y, t, 2);
                            float weight = 1.0f/(completeMatch(x, y, t, 3) + 1);

                            if (sourceMask.defined()) { weight *= sourceMask(x, y, t, 0); }

                            for (int dy = -patchSize/2; dy <= patchSize/2; dy++) {
                                if (y+dy < 0) continue;
                                if (y+dy >= source.height) break;
                                for (int dx = -patchSize/2; dx <= patchSize/2; dx++) {
                                    if (x+dx < 0) continue;
                                    if (x+dx >= source.width) break;

                                    float w = weight;
                                    if (targetMask.defined()) {
                                        w *= targetMask(dstX + dx, dstY + dy, dstT, 0);
                                    }
                                    if (w == 0) continue;

                                    for (int c = 0; c < source.channels; c++) {
                                        out(dstX+dx, dstY+dy, dstT, c) += w*source(x+dx, y+dy, t, c);
                                    }
                                    out(dstX+dx, dstY+dy, dstT, source.channels) += w;
                                }
                            }
                        }
                    }
                }
            }
        }

        if (alpha != 1) {
            // COHERENCE TERM
            Image coherentMatch = PatchMatch::apply(target, source, sourceMask,
                                                    numIterPM, patchSize);
            // For every patch in the target, pull from the nearest match in the source
            for (int t = 0; t < target.frames; t++) {
                for (int y = 0; y < target.height; y++) {
                    for (int x = 0; x < target.width; x++) {

                        float patchWeight = targetWeight.defined() ? targetWeight(x, y, t, 0) : 1;

                        // Target is completely defined here, this iteration is useless
                        if (patchWeight < 1e-10) continue;

                        int dstX = (int)coherentMatch(x, y, t, 0);
                        int dstY = (int)coherentMatch(x, y, t, 1);
                        int dstT = (int)coherentMatch(x, y, t, 2);
                        float weight = 1/(coherentMatch(x, y, t, 3)+1);

                        // Make mostly-defined patches up to 100x more forceful
                        weight *= (1.01 - patchWeight);

                        for (int dy = -patchSize/2; dy <= patchSize/2; dy++) {
                            if (y+dy < 0) { continue; }
                            if (y+dy >= out.height) { break; }
                            for (int dx = -patchSize/2; dx <= patchSize/2; dx++) {
                                if (x+dx < 0) continue;
                                if (x+dx >= out.width) break;
                                float w = weight;
                                if (targetMask.defined()) w *= targetMask(x+dx, y+dy, t, 0);
                                for (int c = 0; c < source.channels; c++) {
                                    out(x+dx, y+dy, t, c) += w*source(dstX+dx, dstY+dy, dstT, c);
                                }
                                out(x+dx, y+dy, t, source.channels) += w;
                            }
                        }
                    }
                }
            }
        }

        // rewrite the target using the homogeneous output

        for (int c = 0; c < out.channels-1; c++) {
            out.channel(c) /= out.channel(out.channels-1) + 1e-10;
        }

        if (targetMask.defined()) {
            Composite::apply(target, out.selectChannels(0, target.channels), targetMask);
        } else {
            target.set(out.selectChannels(0, target.channels));
        }

    }
    printf("\n");
}

void Heal::help() {
    printf("-heal takes an image and a mask, and reconstructs the portion of"
           " the image where the mask is zero using patches from the rest of the"
           " image. It uses the patchmatch algorithm for acceleration. The"
           " arguments include the number of iterations to run per scale, and the"
           " number of iterations of patchmatch to run. Both default to five.\n"
           "\n"
           "Usage: ImageStack -load mask.png -load image.jpg -heal -display\n");
}

bool Heal::test() {
    // Largely the same test as inpaint

    // A circular mask
    Image mask(99, 97, 1, 1);
    Expr::X x; Expr::Y y;
    mask.set(Select((x-50)*(x-50) + (y-50)*(y-50) < 30*30, 0, 1));

    // Patch a corrupted region of a smooth ramp
    Image im(99, 97, 1, 3);
    im.set((x + y)/100);
    Image corrupted = im.copy();
    Noise::apply(corrupted.region(40, 40, 0, 0, 20, 20, 1, 3), -20, 20);
    Heal::apply(corrupted, mask);
    return nearlyEqual(im, corrupted);
}

void Heal::parse(vector<string> args) {
    int numIter = 5;
    int numIterPM = 5;

    assert(args.size() < 3, "-heal takes zero, one, or two arguments\n");

    if (args.size() > 0) { numIter = readInt(args[0]); }
    if (args.size() > 1) { numIterPM = readInt(args[1]); }

    apply(stack(0), stack(1), numIter, numIterPM);
}

void Heal::apply(Image image, Image mask, int numIter, int numIterPM) {

    // First inpaint across the hole and add noise
    image.set(Inpaint::apply(image, mask));
    Image noise(image.width, image.height, image.frames, image.channels);
    Noise::apply(noise, -0.3, 0.3);
    for (int c = 0; c < image.channels; c++) {
        noise.channel(c) *= (1-mask);
    }
    image += noise;

    BidirectionalSimilarity::apply(image.copy(), image, mask, 1-mask, 0, numIter, numIterPM);

}
#include "footer.h"
