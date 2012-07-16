#include "main.h"
#include "Complex.h"
#include "Statistics.h"
#include "header.h"

void ComplexMultiply::help() {
    pprintf("-complexmultiply multiplies the top image in the stack by the second"
            " image in the stack, using 2 \"complex\" images as its input - a"
            " \"complex\" image is one where channel 2*n is the real part of the nth"
            " channel and channel 2*n + 1 is the imaginary part of the nth"
            " channel. Using zero arguments results in a straight multiplication"
            " (a + bi) * (c + di), using one argument results in a conjugate"
            " multiplication (a - bi) * (c + di).\n"
            "\n"
            "Usage: ImageStack -load a.tga -load b.tga -complexmultiply -save out.tga.\n");
}

bool ComplexMultiply::test() {
    Image a(324, 243, 4, 4);
    Image b(324, 243, 4, 4);
    Noise::apply(a, -1, 1);

    // a * conj(b) == b * conj(a)
    Image c = a.copy();
    ComplexMultiply::apply(c, b, true);
    Image d = b.copy();
    ComplexMultiply::apply(d, a, true);
    if (!nearlyEqual(c, d)) return false;

    return true;
}

void ComplexMultiply::parse(vector<string> args) {
    assert(args.size() < 2, "-complexmultiply takes zero or one arguments\n");
    if (stack(0).channels == 2 && stack(1).channels > 2) {
        apply(stack(1), stack(0), (bool)args.size());
        pop();
    } else {
        apply(stack(0), stack(1), (bool)args.size());
        pull(1);
        pop();
    }
}

void ComplexMultiply::apply(Image a, Image b, bool conj) {
    assert(a.channels % 2 == 0 && b.channels % 2 == 0,
           "-complexmultiply requires images with an even number of channels (%d %d)\n",
           a.channels, b.channels);

    assert(a.frames == b.frames &&
           a.width == b.width &&
           a.height == b.height,
           "images must be the same size\n");

    if (a.channels == 2 && b.channels == 2) {
        // Scalar times scalar
        Image a_real = a.channel(0), a_imag = a.channel(1);
        Image b_real = b.channel(0), b_imag = b.channel(1);
        if (conj) {
            a.setChannels(a_real*b_real + a_imag*b_imag,
                          a_imag*b_real - a_real*b_imag);
        } else {
            a.setChannels(a_real*b_real - a_imag*b_imag,
                          a_imag*b_real + a_real*b_imag);
        }
    } else if (b.channels == 2) {
        // Vector times scalar
        for (int c = 0; c < a.channels; c += 2) {
            apply(a.selectChannels(c, 2), b, conj);
        }
    } else {
        // Vector times vector (elementwise)
        for (int c = 0; c < a.channels; c += 2) {
            apply(a.selectChannels(c, 2),
                  b.selectChannels(c, 2), conj);
        }
    }
}


void ComplexDivide::help() {
    pprintf("-complexdivide divides the top image in the stack by the second image"
            " in the stack, using 2 \"complex\" images as its input - a \"complex\""
            " image is one where channel 2*n is the real part of the nth channel and"
            " channel 2*n + 1 is the imaginary part of the nth channel. Using zero"
            " arguments results in a straight division (a + bi) / (c + di). Using"
            " one argument results in a conjugate division (a - bi) / (c + di).\n"
            "\n"
            "Usage: ImageStack -load a.tga -load b.tga -complexdivide -save out.tga.\n");
}

bool ComplexDivide::test() {
    Image a(123, 234, 4, 2);
    Image b(123, 234, 4, 2);
    Image c(123, 234, 4, 2);

    // (a + b) / (conj(conj(c))) = a / c + b / c
    Noise::apply(a, -1, 1);
    Noise::apply(b, -1, 1);
    Noise::apply(c, 1, 2);
    Image d = a + b;
    Image cc = c.copy();
    ComplexConjugate::apply(cc);
    ComplexDivide::apply(d, cc, true);
    ComplexDivide::apply(a, c);
    ComplexDivide::apply(b, c);
    return nearlyEqual(d, a+b);
}

void ComplexDivide::parse(vector<string> args) {
    assert(args.size() == 0 || args.size() == 1,
           "-complexdivide takes zero or one arguments\n");
    if (stack(0).channels == 2 && stack(1).channels > 2) {
        apply(stack(1), stack(0), (bool)args.size());
        pop();
    } else {
        apply(stack(0), stack(1), (bool)args.size());
        pull(1);
        pop();
    }
}

void ComplexDivide::apply(Image a, Image b, bool conj) {
    assert(a.channels % 2 == 0 && b.channels % 2 == 0,
           "-complexdivide requires images with an even number of channels\n");

    assert(a.frames == b.frames &&
           a.width == b.width &&
           a.height == b.height,
           "images must be the same size\n");

    if (a.channels == 2 && b.channels == 2) {
        // Scalar over scalar
        Image a_real = a.channel(0), a_imag = a.channel(1);
        Image b_real = b.channel(0), b_imag = b.channel(1);
        auto denom = b_real*b_real + b_imag*b_imag;
        if (conj) {
            a.setChannels((a_real * b_real - a_imag * b_imag) / denom,
                          (a_imag * b_real + a_real * b_imag) / denom);
        } else {
            a.setChannels((a_real * b_real + a_imag * b_imag) / denom,
                          (a_imag * b_real - a_real * b_imag) / denom);
        }
    } else if (b.channels == 2) {
        // Vector over scalar
        for (int c = 0; c < a.channels; c += 2) {
            apply(a.selectChannels(c, 2), b, conj);
        }
    } else {
        // Vector over vector (elementwise)
        for (int c = 0; c < a.channels; c += 2) {
            apply(a.selectChannels(c, 2),
                  b.selectChannels(c, 2), conj);
        }
    }
}


void ComplexReal::help() {
    pprintf("-complexreal takes a \"complex\" image, in which the even channels"
            " represent the real component and the odd channels represent the"
            " imaginary component, and produces an image containing only the real"
            " channels.\n"
            "\n"
            "Usage: ImageStack -load a.png -fftreal -complexreal -display\n");
}

bool ComplexReal::test() {
    // a * conj(a) == re(a)^2 + imag(a)^2 == |a|^2
    Image a(123, 234, 4, 2);
    Noise::apply(a, -1, 1);
    Image a_real = ComplexReal::apply(a);
    Image a_imag = ComplexImag::apply(a);
    Image mag = ComplexMagnitude::apply(a);
    ComplexMultiply::apply(a, a, true);
    
    if (!nearlyEqual(ComplexReal::apply(a), a_real*a_real + a_imag*a_imag)) return false;
    if (!nearlyEqual(ComplexReal::apply(a), mag*mag)) return false;

    return true;
}

void ComplexReal::parse(vector<string> args) {
    assert(args.size() == 0, "-complexreal takes no arguments\n");
    Image im = apply(stack(0));
    pop();
    push(im);
}

Image ComplexReal::apply(Image im) {
    assert(im.channels % 2 == 0,
           "complex images must have an even number of channels\n");

    Image out(im.width, im.height, im.frames, im.channels/2);

    for (int c = 0; c < out.channels; c++) {
        out.channel(c).set(im.channel(2*c));
    }

    return out;
}

void RealComplex::help() {
    pprintf("-realcomplex takes a \"real\" image, and converts it to a \"complex\""
            " image, in which the even channels represent the real component and"
            " the odd channels represent the imaginary component.\n"
            "\n"
            "Usage: ImageStack -load a.png -realcomplex -fft -display\n");
}

bool RealComplex::test() {
    Image a(123, 234, 3, 2);
    Noise::apply(a, -1, 1);
    Image b = RealComplex::apply(a);
    Image c = ComplexReal::apply(b);
    Stats s(ComplexImag::apply(b));
    return (nearlyEqual(c, a) && s.mean() == 0 && s.variance() == 0);
}

void RealComplex::parse(vector<string> args) {
    assert(args.size() == 0, "-complexreal takes no arguments\n");
    Image im = apply(stack(0));
    pop();
    push(im);
}

Image RealComplex::apply(Image im) {
    Image out(im.width, im.height, im.frames, im.channels*2);

    for (int c = 0; c < im.channels; c++) {
        out.channel(2*c).set(im.channel(c));
    }

    return out;
}

void ComplexImag::help() {
    pprintf("-compleximag takes a \"complex\" image, in which the even channels"
            " represent the real component and the odd channels represent the"
            " imaginary component, and produces an image containing only the imaginary"
            " channels.\n"
            "\n"
            "Usage: ImageStack -load a.png -fftreal -compleximag -display\n");
}

bool ComplexImag::test() {
    // tested by ComplexReal
    return true;
}

void ComplexImag::parse(vector<string> args) {
    assert(args.size() == 0, "-compleximag takes no arguments\n");
    Image im = apply(stack(0));
    pop();
    push(im);
}

Image ComplexImag::apply(Image im) {
    assert(im.channels % 2 == 0,
           "complex images must have an even number of channels\n");

    Image out(im.width, im.height, im.frames, im.channels/2);

    for (int c = 0; c < out.channels; c++) {
        out.channel(c).set(im.channel(2*c+1));
    }

    return out;
}


void ComplexMagnitude::help() {
    pprintf("-complexmagnitude takes a \"complex\" image, in which the even channels"
            " represent the real component and the odd channels represent the"
            " imaginary component, and produces an image containing the complex"
            " magnitude\n"
            "\n"
            "Usage: ImageStack -load a.png -fftreal -complexmagnitude -display\n");
}

bool ComplexMagnitude::test() {
    // tested by ComplexReal
    return true;
}

void ComplexMagnitude::parse(vector<string> args) {
    assert(args.size() == 0, "-complexmagnitude takes no arguments\n");
    Image im = apply(stack(0));
    pop();
    push(im);
}

Image ComplexMagnitude::apply(Image im) {
    assert(im.channels % 2 == 0,
           "complex images must have an even number of channels\n");

    Image out(im.width, im.height, im.frames, im.channels/2);

    for (int c = 0; c < out.channels; c++) {
        Image real = im.channel(2*c);
        Image imag = im.channel(2*c+1);
        out.channel(c).set(sqrt(real*real + imag*imag));
    }

    return out;
}



void ComplexPhase::help() {
    pprintf("-complexphase takes a \"complex\" image, in which the even channels"
            " represent the real component and the odd channels represent the"
            " imaginary component, and produces an image containing the complex"
            " phase\n"
            "\n"
            "Usage: ImageStack -load a.png -fftreal -complexphase -display\n");
}

bool ComplexPhase::test() {
    Image a(123, 234, 3, 2);
    Noise::apply(a, 1, 2);
    a = RealComplex::apply(a);

    // Multiply real-valued image by 1+i and check the phase
    Image b(123, 234, 3, 2);
    b.set(1);
    ComplexMultiply::apply(a, b);
    b = ComplexPhase::apply(a);
    Stats s(b);
    if (!(nearlyEqual(s.mean(), M_PI/4) &&
          nearlyEqual(s.variance(), 0))) return false;

    // Squaring should double phase
    a.set(0);
    Noise::apply(a, 1, 2);
    b = a.copy();
    ComplexMultiply::apply(b, b);
    return nearlyEqual(ComplexPhase::apply(b),
                       2*ComplexPhase::apply(a));
}

void ComplexPhase::parse(vector<string> args) {
    assert(args.size() == 0, "-complexphase takes no arguments\n");
    Image im = apply(stack(0));
    pop();
    push(im);
}

Image ComplexPhase::apply(Image im) {
    assert(im.channels % 2 == 0, "complex images must have an even number of channels\n");

    Image out(im.width, im.height, im.frames, im.channels/2);

    for (int c = 0; c < out.channels; c++) {
        Image real = im.channel(2*c);
        Image imag = im.channel(2*c+1);
        out.channel(c).set(Expr::atan2(imag, real));
    }

    return out;
}


void ComplexConjugate::help() {
    pprintf("-complexconjugate takes a \"complex\" image, in which the even channels"
            " represent the real component and the odd channels represent the"
            " imaginary component, and produces an image containing the complex"
            " conjugate\n"
            "\n"
            "Usage: ImageStack -load a.png -fftreal -complexconjugate -display\n");
}

void ComplexConjugate::parse(vector<string> args) {
    assert(args.size() == 0, "-complexconjugate takes no arguments\n");
    apply(stack(0));
}

bool ComplexConjugate::test() {
    // Tested by complex multiply
    return true;
}

void ComplexConjugate::apply(Image im) {
    assert(im.channels % 2 == 0, "complex images must have an even number of channels\n");

    for (int c = 1; c < im.channels; c+=2) {
        im.channel(c).set(-im.channel(c));
    }
}

#include "footer.h"
