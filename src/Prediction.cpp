#include "main.h"
#include "Prediction.h"
#include "Filter.h"
#include "Paint.h"
#include "File.h"
#include "Calculus.h"
#include "Statistics.h"
#include "Convolve.h"
#include "header.h"

void Inpaint::help() {
    printf("\n-inpaint takes the image on the top of the stack, and a one channel mask of the\n"
           "same size second on the stack, and diffuses areas of the image where the mask is\n"
           "high into areas of the image where the mask is low. Image pixels with mask of 1\n"
           "are unchanged.\n\n"
           "Usage: ImageStack -push 1 640 480 1 -eval \"(X > 0.5)*(X < 0.6)\" -load in.jpg\n"
           "                  -inpaint -save out.jpg\n\n");
}

bool Inpaint::test() {
    // A circular mask
    Image mask(99, 97, 1, 1);
    Expr::X x; Expr::Y y;
    mask.set(Select((x-50)*(x-50) + (y-50)*(y-50) < 30*30, 0, 1));

    // Patch a corrupted region of a smooth ramp
    Image im(99, 97, 1, 3);
    im.set((x + y)/100.0f);
    Image corrupted = im.copy();
    Noise::apply(corrupted.region(40, 40, 0, 0, 20, 20, 1, 3), -20, 20);
    Image after = Inpaint::apply(corrupted, mask);
    if (!nearlyEqual(im, after)) return false;

    // There should be no strong gradients in the output
    GradMag::apply(after);
    after = after.region(10, 10, 0, 0, 80, 80, 1, 3);
    Stats s(after);
    if (s.maximum() > 0.01) return false;

    // Smooth a hole within a noise image
    im.set(0);
    Noise::apply(im, 0, 1);
    after = Inpaint::apply(im, mask);

    // Statistics within the hole should be mean of 0.5, zero variance
    s = Stats(after.region(45, 45, 0, 0, 10, 10, 1, 3));
    if (s.mean() < 0.45 ||
        s.mean() > 0.55 ||
        !nearlyEqual(s.variance(), 0)) return false;

    // Outside the hole should be untouched
    for (int c = 0; c < 3; c++) {
        im.channel(c) *= mask;
        after.channel(c) *= mask;
    }
    return nearlyEqual(im, after);
}

void Inpaint::parse(vector<string> args) {
    assert(args.size() == 0, "-inpaint takes no arguments\n");
    Image im = apply(stack(0), stack(1));
    pop();
    push(im);
}

Image Inpaint::apply(Image im, Image mask) {
    assert(im.width == mask.width &&
           im.height == mask.height &&
           im.frames == mask.frames,
           "mask must be the same size as the image\n");
    assert(mask.channels == 1,
           "mask must have one channel\n");

    // Make a mask that just selects the outline
    Image boundaryMask = mask.copy();
    FastBlur::apply(boundaryMask, 1, 1, 1);
    boundaryMask.set(max(mask - boundaryMask, 0));

    // Mask out the input image with it
    Image boundary = im.copy();
    for (int c = 0; c < im.channels; c++) {
        boundary.channel(c) *= boundaryMask;
    }

    Image blurred = boundary.copy();
    Image blurredMask = boundaryMask.copy();

    const int J = 10;
    for (int i = J; i >= 0; i--) {
        float alpha = powf(i, 2)/4;
        // Smooth
        FastBlur::apply(blurred, alpha, alpha, alpha);
        FastBlur::apply(blurredMask, alpha, alpha, alpha);
        // Reimpose boundary conditions
        Composite::apply(blurred, boundary, boundaryMask);
        Composite::apply(blurredMask, boundaryMask, boundaryMask);
    }

    // Normalize
    for (int c = 0; c < im.channels; c++) {
        blurred.channel(c) /= blurredMask;
    }

    // Reattach the rest of the image
    Composite::apply(blurred, im, mask);
    return blurred;
}



void SeamlessClone::help() {
    pprintf("-seamlessclone composites the top image in the stack over the next image in"
            " the stack, using the last channel in the top image in the stack as alpha."
            " The composite is done in such a way as to avoid hard edges around the"
            " boundaries. If the top image in the stack has only one channel, it"
            " interprets this as a mask, and composites the second image in the"
            " stack over the third image in the stack using that mask.\n"
            "\n"
            "Usage: ImageStack -load a.jpg -load b.jpg -load mask.png -seamlessclone\n"
            "       ImageStack -load a.jpg -load b.jpg -evalchannels [0] [1] [2] \\\n"
            "       \"x>width/2\" -seamlessclone -display\n\n");
}

bool SeamlessClone::test() {
    // This op is a trivial use of inpaint. If inpaint works so does
    // this. Let's just make sure it doesn't crash.


    // A circular mask
    Image mask(99, 97, 1, 1);
    Expr::X x; Expr::Y y;
    mask.set(Select((x-50)*(x-50) + (y-50)*(y-50) < 30*30, 0, 1));

    // The background is a smooth ramp
    Image background(99, 97, 1, 3);
    background.set((x + 2*y)/300);

    // The foreground to paste on is noise
    Image foreground(99, 97, 1, 3);
    Noise::apply(foreground, 0, 1);

    Image im = background.copy();
    SeamlessClone::apply(im, foreground, mask);

    // The laplacian of the result should be as if you just pasted the
    // laplacian of the foreground on top of the background.
    Image lap(3, 3, 1, 1);
    lap(1, 0) = lap(0, 1) = lap(2, 1) = lap(1, 2) = 1;
    lap(1, 1) = -4;
    Image lapFg = Convolve::apply(foreground, lap, Convolve::Clamp);
    Image lapBg = Convolve::apply(background, lap, Convolve::Clamp);
    Image lapIm = Convolve::apply(im, lap, Convolve::Clamp);
    Composite::apply(lapBg, lapFg, 1-mask);
    return nearlyEqual(lapIm, lapBg);
}

void SeamlessClone::parse(vector<string> args) {
    assert(args.size() == 0, "-seamlessclone takes no arguments\n");

    if (stack(0).channels == 1) {
        apply(stack(2), stack(1), stack(0));
        pop();
        pop();
    } else {
        apply(stack(1), stack(0));
        pop();
    }
}

void SeamlessClone::apply(Image dst, Image src) {
    assert(src.channels > 1, "Source image needs at least two channels\n");
    assert(src.channels == dst.channels || src.channels == dst.channels + 1,
           "Source image and destination image must either have matching channel"
           " counts (if they both have an alpha channel), or the source image"
           " should have one more channel than the destination.\n");
    assert(dst.frames == src.frames && dst.width == src.width
           && dst.height == src.height,
           "The source and destination images must be the same size\n");

    if (src.channels > dst.channels) {
        apply(dst,
              src.region(0, 0, 0, 0,
                         src.width, src.height,
                         src.frames, dst.channels),
              src.channel(dst.channels));

    } else {
        apply(dst, src, src.channel(dst.channels-1));
    }
}

void SeamlessClone::apply(Image dst, Image src, Image mask) {
    assert(src.channels == dst.channels, "The source and destination images must have the same number of channels\n");

    assert(dst.frames == src.frames && dst.width == src.width && dst.height == src.height,
           "The source and destination images must be the same size\n");
    assert(dst.frames == mask.frames && dst.width == mask.width && dst.height == mask.height,
           "The source and destination images must be the same size as the mask\n");

    assert(mask.channels == 1, "Mask must have one channel\n");

    // Generate a smooth patch to fix discontinuities between source and destination
    Image patch = Inpaint::apply(dst-src, mask);

    for (int c = 0; c < dst.channels; c++) {
        dst.channel(c).set((1-mask)*(src.channel(c) + patch.channel(c)) + mask*dst.channel(c));
    }
}

#include "footer.h"
