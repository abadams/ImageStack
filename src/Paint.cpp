#include "main.h"
#include "Paint.h"
#include "Parser.h"
#include "File.h"
#include "Geometry.h"
#include "header.h"

using namespace Expr;

void Eval::help() {
    printf("\n-eval takes a simple expression and evaluates it, writing the result to the\n"
           "current image.\n\n");
    Expression::help();
    printf("Usage: ImageStack -push 128 128 128 1 -eval \"(x*y*t)^0.5\" -save out.tga\n\n");
}

bool Eval::test() {
    // This also tests much of the Expr module

    Image im(123, 234, 2, 3);
    Noise::apply(im, 0, 1);

    // Test basic arithmetic and channel selection
    {
        printf("Testing arithmetic\n");
        Image a = Eval::apply(im, "[1]*2 + val");
        Image b(a.width, a.height, a.frames, a.channels);
        b.setChannels(im.channel(1)*2 + im.channel(0),
                      im.channel(1)*2 + im.channel(1),
                      im.channel(1)*2 + im.channel(2));
        if (!nearlyEqual(a, b)) return false;
    }

    // Test some transcendentals
    {
        printf("Testing transcendentals\n");
        Image a = Eval::apply(im, "abs(cos(50*(log(exp(val-2) + 1))))");
        Image b = abs(cos(50*(log(exp(im-2) + 1))));
        if (!nearlyEqual(a, b)) return false;
    }


    // Test conditional
    {
        printf("Testing conditionals\n");
        Stats ims(im);
        Image a = Eval::apply(im, "val > 0.5 ? (val / mean(2) + skew(0)) : -covariance(0, 1)");
        Image b = Select(im > 0.5, im / ims.mean(2) + ims.skew(0), -ims.covariance(0, 1));
        if (!nearlyEqual(a, b)) return false;
    }

    // Test some more comparisons
    {
        printf("Testing comparisons\n");
        Stats ims(im);
        Image a = Eval::apply(im, "(val > 0.5) + (x <= 50) + (y >= 10) + (c < 2) + (y != 23) + (x == 22)");
        Image b = (Select(im > 0.5, 1.0f, 0.0f) +
                   Select(X() <= 50, 1.0f, 0.0f) +
                   Select(Y() >= 10, 1.0f, 0.0f) +
                   Select(C() < 2, 1.0f, 0.0f) +
                   Select(Y() != 23, 1.0f, 0.0f) +
                   Select(X() == 22, 1.0f, 0.0f));
        if (!nearlyEqual(a, b)) return false;
    }

    // Test sampling
    {
        printf("Testing sampling\n");
        Stats ims(im);
        Image a = Eval::apply(im, "[x*0.8 + y*0.2, -x*0.2 + y*0.8]");
        float matrix[] = {0.8, 0.2, 0, -0.2, 0.8, 0};
        Image b = AffineWarp::apply(im, matrix);
        if (!nearlyEqual(a, b)) return false;
    }

    return true;
}

void Eval::parse(vector<string> args) {
    assert(args.size() == 1, "-eval takes exactly one argument\n");
    Image im = apply(stack(0), args[0]);
    pop();
    push(im);
}

Image Eval::apply(Image im, string expression_) {
    Expression expression(expression_);
    Expression::State state(im);

    Image out(im.width, im.height, im.frames, im.channels);

    for (state.c = 0; state.c < im.channels; state.c++) {
        for (state.t = 0; state.t < im.frames; state.t++) {
            for (state.y = 0; state.y < im.height; state.y++) {
                for (state.x = 0; state.x < im.width; state.x++) {
                    out(state.x, state.y, state.t, state.c) = expression.eval(state);
                }
            }
        }
    }

    return out;
}



void EvalChannels::help() {
    printf("\n-evalchannels takes some expressions and evaluates them, writing the results\n"
           "to an image with that many channels.\n\n");
    Expression::help();
    printf("Usage: ImageStack -push 128 128 128 1 -evalchannels \"[0]*2\" \"[1]*2 + [0]\"\n"
           "                  -save out.tga\n\n");
}

bool EvalChannels::test() {
    // Largely already tested by -eval
    vector<string> expressions(3);
    expressions[0] = "[2] + 1";
    expressions[1] = "[0] / 17";
    expressions[2] = "[1] + (x > 10 ? 50 : 0)";
    Image im(123, 234, 2, 3);
    Noise::apply(im, 0, 1);
    Image a = EvalChannels::apply(im, expressions);
    Image b = im.copy();
    b.setChannels(im.channel(2)+1,
                  im.channel(0)/17,
                  im.channel(1) + Select(X() > 10, 50, 0));
    return nearlyEqual(a, b);
}

void EvalChannels::parse(vector<string> args) {
    Image im = apply(stack(0), args);
    pop();
    push(im);
}

Image EvalChannels::apply(Image im, vector<string> expressions_) {
    vector<Expression *> expressions(expressions_.size());
    for (size_t i = 0; i < expressions_.size(); i++) {
        expressions[i] = new Expression(expressions_[i]);
    }

    int channels = (int)expressions_.size();

    Image out(im.width, im.height, im.frames, channels);

    Expression::State state(im);

    for (state.c = 0; state.c < channels; state.c++) {
        for (state.t = 0; state.t < im.frames; state.t++) {
            for (state.y = 0; state.y < im.height; state.y++) {
                for (state.x = 0; state.x < im.width; state.x++) {
                    out(state.x, state.y, state.t, state.c) = expressions[state.c]->eval(state);
                }
            }
        }
    }

    for (size_t i = 0; i < expressions.size(); i++) { delete expressions[i]; }

    return out;
}

void Plot::help() {
    printf("\n-plot takes images with height 1 and range [0, 1], and graphs them.\n"
           "It takes three arguments: the width and height of the resulting graph,\n"
           "and the line thickness to use for the plot. The resulting graph will\n"
           "have the same number of frames and channels as the input.\n\n");
}

bool Plot::test() {
    Image a(123, 1, 1, 7);
    Image p = Plot::apply(a, 512, 512, 3);
    if (p.channels != 7 ||
        p.width != 512 ||
        p.height != 512) return false;
    Stats s(p);
    return (s.mean() > 0.0001 && s.mean() < 0.02 && s.variance() < 0.001);
}

void Plot::parse(vector<string> args) {
    Image im = apply(stack(0), readInt(args[0]), readInt(args[1]), readFloat(args[2]));
    push(im);
}

Image Plot::apply(Image im, int width, int height, float lineThickness) {
    Image out(width, height, im.frames, im.channels);

    // convert from diameter to radius
    lineThickness /= 2;

    float widthScale = (float)out.width / im.width;

    for (int t = 0; t < im.frames; t++) {
        for (int i = 0; i < im.width-1; i++) {
            for (int c = 0; c < im.channels; c++) {
                float x1 = i*widthScale;
                float x2 = (i+1)*widthScale;
                float y1 = ((1-im(i, 0, t, c)) * out.height + 0.5);
                float y2 = ((1-im(i+1, 0, t, c)) * out.height + 0.5);
                int minY, maxY;
                int minX = (int)floorf(x1 - lineThickness - 1);
                int maxX = (int)ceilf(x2 + lineThickness + 1);

                if (y1 < y2) {
                    minY = (int)floorf(y1 - lineThickness - 1);
                    maxY = (int)ceilf(y2 + lineThickness + 1);
                } else {
                    minY = (int)floorf(y2 - lineThickness - 1);
                    maxY = (int)ceilf(y1 + lineThickness + 1);
                }

                float deltaX = x2 - x1;
                float deltaY = y2 - y1;
                float segmentLength = sqrtf(deltaX * deltaX + deltaY * deltaY);
                deltaX /= segmentLength;
                deltaY /= segmentLength;

                for (int y = minY; y <= maxY; y++) {
                    if (y < 0 || y >= out.height) { continue; }
                    for (int x = minX; x <= maxX; x++) {
                        if (x < 0 || x >= out.width) { continue; }
                        float bestDistance = lineThickness+2;

                        // check distance to the points
                        float d = (x1 - x)*(x1 - x) + (y1 - y)*(y1 - y);
                        if (d < bestDistance*bestDistance) { bestDistance = sqrtf(d); }
                        if (i == im.width-2) {
                            // check the last point
                            d = sqrtf((x2 - x)*(x2 - x) + (y2 - y)*(y2 - y));
                            if (d < bestDistance) { bestDistance = d; }
                        }

                        // check distance to the line
                        float alpha = deltaX*(x-x1) + deltaY*(y-y1);
                        float beta  = -deltaY*(x-x1) + deltaX*(y-y1);
                        if (alpha > 0 && alpha < segmentLength) {
                            if (beta < 0) { beta = -beta; }
                            if (beta < bestDistance) { bestDistance = beta; }
                        }

                        float result = 0;
                        if (bestDistance < lineThickness - 0.5) { result = 1; }
                        else if (bestDistance < lineThickness + 0.5) { result = (lineThickness + 0.5) - bestDistance; }
                        if (out(x, y, t, c) < result) { out(x, y, t, c) = result; }
                    }
                }
            }
        }
    }

    return out;
}



void Composite::help() {
    pprintf("-composite composites the top image in the stack over the next image"
            " in the stack, using the last channel in the top image in the stack as"
            " alpha. If the top image in the stack has only one channel, it"
            " interprets this as a mask, and composites the second image in the"
            " stack over the third image in the stack using that mask.\n"
            "\n"
            "Usage: ImageStack -load a.jpg -load b.jpg -load mask.png -composite\n"
            "       ImageStack -load a.jpg -load b.jpg -evalchannels [0] [1] [2] \\\n"
            "       \"x>width/2\" -composite -display\n\n");
}

bool Composite::test() {
    Image mask(123, 234, 1, 1);

    mask.set((X()+Y())/(123.0f+234.0f));

    Image a(123, 234, 1, 3);
    Image b(123, 234, 1, 3);
    Noise::apply(a, 0, 6);
    Noise::apply(b, -3, 2);
    Image c = a.copy();
    Composite::apply(c, b, mask);

    for (int i = 0; i < 100; i++) {
        int y = randomInt(0, a.height-1);
        int x = randomInt(0, a.width-1);
        float m = (x+y)/(123.0f + 234.0f);
        float val = c(x, y, 1);
        float correct = m*b(x, y, 1) + (1-m)*a(x, y, 1);
        if (!nearlyEqual(val, correct)) return false;
    }

    return true;
}

void Composite::parse(vector<string> args) {
    assert(args.size() == 0, "-composite takes no arguments\n");

    if (stack(0).channels == 1) {
        apply(stack(2), stack(1), stack(0));
        pop();
        pop();
    } else {
        apply(stack(1), stack(0));
        pop();
    }
}

void Composite::apply(Image dst, Image src) {
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

void Composite::apply(Image dst, Image src, Image mask) {
    assert(src.channels == dst.channels, "The source and destination images must have the same number of channels\n");

    assert(dst.frames == src.frames && dst.width == src.width && dst.height == src.height,
           "The source and destination images must be the same size\n");
    assert(dst.frames == mask.frames && dst.width == mask.width && dst.height == mask.height,
           "The source and destination images must be the same size as the mask\n");

    for (int c = 0; c < dst.channels; c++) {
        dst.channel(c).set(mask*src.channel(c) + (1-mask)*dst.channel(c));
    }
}

#include "footer.h"
