#ifndef NO_FFTW
#include "main.h"
#include "DFT.h"
#include "Geometry.h"
#include "Stack.h"
#include "Arithmetic.h"
#include "Complex.h"
#include "Display.h"
#include "Statistics.h"
#include "Calculus.h"
#include "File.h"
#include <fftw3.h>
#include "header.h"

void DCT::help() {
    pprintf("-dct performs a real discrete cosine transform on the current"
            " image, over the dimensions given in the argument. The signal is"
            " treated as having reflective symmetry. If no arguments are"
            " given, every dimension is transformed.\n"
            "\n"
            "Usage: ImageStack -load a.png -dct xy -save freq.png\n");

}

bool DCT::test() {
    Image a(123, 234, 10, 3);
    Noise::apply(a, 0, 1);
    Image b = a.copy();
    // Check it's a self-inverse
    DCT::apply(a, false, true, true);
    DCT::apply(a, true, false, true);
    DCT::apply(a, true, true, false);
    if (!nearlyEqual(b, a)) return false;

    // Check it's an orthogonal transform
    b.set(0); Noise::apply(b, 0, 1);
    double d1 = Stats(a - b).variance();
    DCT::apply(a, true, true, false);
    DCT::apply(b, true, true, false);
    double d2 = Stats(a - b).variance();
    if (!nearlyEqual(d1, d2)) return false;

    // Check a single cosine curve creates a spike
    a.set(cos(16 * M_PI * Expr::X() / 122.0));
    DCT::apply(a, true, true, true);
    if (a(16, 0, 0, 0) < 1) return false;
    a(16, 0, 0, 0) = 0;
    a(16, 0, 0, 1) = 0;
    a(16, 0, 0, 2) = 0;
    a *= a;
    return Stats(a).maximum() < 0.0001;
}

void DCT::parse(vector<string> args) {
    assert(args.size() < 2, "-dct takes zero or one argument\n");

    bool x = true, y = true, t = true;
    if (args.size() == 1) {
        x = y = t = false;
        for (size_t i = 0; i < args[0].size(); i++) {
            switch (args[0][i]) {
            case 'x':
                x = true;
                break;
            case 'y':
                y = true;
                break;
            case 't':
                t = true;
                break;
            default:
                panic("Unknown dimension: %c\n", args[0][i]);
                break;
            }
        }

    }

    apply(stack(0), x, y, t);
}

void DCT::apply(Image im, bool transformX, bool transformY, bool transformT) {
    if (im.width == 1) { transformX = false; }
    if (im.height == 1) { transformY = false; }
    if (im.frames == 1) { transformT = false; }

    // rank 0
    if (!transformX && !transformY && !transformT) { return; }

    vector<fftwf_iodim> loop_dims;
    vector<fftwf_iodim> fft_dims;

    // X
    {
        fftwf_iodim d = {im.width, 1, 1};
        if (transformX) fft_dims.push_back(d);
        else loop_dims.push_back(d);
    }

    // Y
    {
        fftwf_iodim d = {im.height, im.ystride, im.ystride};
        if (transformY) fft_dims.push_back(d);
        else loop_dims.push_back(d);
    }

    // T
    {
        fftwf_iodim d = {im.frames, im.tstride, im.tstride};
        if (transformT) fft_dims.push_back(d);
        else loop_dims.push_back(d);
    }

    // C
    {
        fftwf_iodim d = {im.channels, im.cstride, im.cstride};
        loop_dims.push_back(d);
    }

    vector<fftw_r2r_kind> kinds(fft_dims.size(), FFTW_REDFT00);

    fftwf_plan plan = fftwf_plan_guru_r2r((int)fft_dims.size(), &fft_dims[0],
                                          (int)loop_dims.size(), &loop_dims[0],
                                          im.baseAddress(), im.baseAddress(),
                                          &kinds[0], FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);

    float m = 1.0;
    if (transformX) m *= 2*(im.width-1);
    if (transformY) m *= 2*(im.height-1);
    if (transformT) m *= 2*(im.frames-1);
    im /= sqrtf(m);
}

void FFT::help() {
    pprintf("-fft performs a fast dft on the current image, whose values"
            " are interpreted as complex. The input is an image with 2*c channels,"
            " where channel 2*i is the real part of the i\'th channel, and channel"
            " 2*i+1 is the imaginary part of the i'th channel. The output image is"
            " laid out the same way.\n"
            "\n"
            "Usage: ImageStack -load a.tmp -fftcomplex -save freq.tmp\n\n");

}


bool FFT::test() {
    Image a(123, 234, 10, 4);
    Noise::apply(a, 0, 1);
    Image b = a.copy();
    // Check it's a self-inverse
    FFT::apply(a, true, true, true);
    IFFT::apply(a, true, true, true);
    if (!nearlyEqual(b, a)) return false;

    // Check it's an orthogonal transform
    ComplexMultiply::apply(b, a, true);
    double d1 = Stats(b).sum();
    FFT::apply(a, true, true, false);
    a /= sqrtf(123*234);
    b = a.copy();
    ComplexMultiply::apply(b, a, true);
    double d2 = Stats(b).sum();
    if (!nearlyEqual(d1, d2)) return false;

    // Check a single shifted curve creates a spike
    a.channel(0).set(cos(16 * M_PI * Expr::X() / 123.0 + M_PI/8));
    a.channel(1).set(sin(16 * M_PI * Expr::X() / 123.0 + M_PI/8));
    a.channel(2).set(a.channel(0));
    a.channel(3).set(a.channel(1));
    FFT::apply(a, true, true, true);
    a /= sqrtf(123*234*10);
    float r = a(8, 0, 0, 0);
    float c = a(8, 0, 0, 1);
    if (r*r + c*c < 1) return false;
    if (!nearlyEqual(atan2(c, r), M_PI/8)) return false;
    a(8, 0, 0, 0) = 0;
    a(8, 0, 0, 1) = 0;
    a(8, 0, 0, 2) = 0;
    a(8, 0, 0, 3) = 0;
    ComplexMultiply::apply(a, a, true);
    return Stats(a).maximum() < 0.0001;
}

void FFT::parse(vector<string> args) {
    assert(args.size() < 2, "-fft takes zero or one argument\n");

    bool x = true, y = true, t = true;
    if (args.size() == 1) {
        x = y = t = false;
        for (size_t i = 0; i < args[0].size(); i++) {
            switch (args[0][i]) {
            case 'x':
                x = true;
                break;
            case 'y':
                y = true;
                break;
            case 't':
                t = true;
                break;
            default:
                panic("Unknown dimension: %c\n", args[0][i]);
                break;
            }
        }

    }

    apply(stack(0), x, y, t);
}

void FFT::apply(Image im, bool transformX, bool transformY, bool transformT, bool inverse) {
    assert(im.channels % 2 == 0, "-fft requires an image with an even number of channels\n");

    if (im.width == 1) { transformX = false; }
    if (im.height == 1) { transformY = false; }
    if (im.frames == 1) { transformT = false; }

    // rank 0
    if (!transformX && !transformY && !transformT) { return; }

    vector<fftwf_iodim> loop_dims;
    vector<fftwf_iodim> fft_dims;

    // X
    {
        fftwf_iodim d = {im.width, 1, 1};
        if (transformX) fft_dims.push_back(d);
        else loop_dims.push_back(d);
    }

    // Y
    {
        fftwf_iodim d = {im.height, im.ystride, im.ystride};
        if (transformY) fft_dims.push_back(d);
        else loop_dims.push_back(d);
    }

    // T
    {
        fftwf_iodim d = {im.frames, im.tstride, im.tstride};
        if (transformT) fft_dims.push_back(d);
        else loop_dims.push_back(d);
    }

    // C
    {
        fftwf_iodim d = {im.channels/2, im.cstride*2, im.cstride*2};
        loop_dims.push_back(d);
    }

    // An inverse fft can be done by swapping real and imaginary parts
    int real_c = inverse ? 1 : 0;
    int imag_c = inverse ? 0 : 1;

    fftwf_plan plan = fftwf_plan_guru_split_dft((int)fft_dims.size(), &fft_dims[0],
                                                (int)loop_dims.size(), &loop_dims[0],
                                                &(im(0, 0, 0, real_c)), &(im(0, 0, 0, imag_c)),
                                                &(im(0, 0, 0, real_c)), &(im(0, 0, 0, imag_c)),
                                                FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);


    if (inverse) {
        float m = 1.0;
        if (transformX) m *= im.width;
        if (transformY) m *= im.height;
        if (transformT) m *= im.frames;
        im /= m;
    }
}


void IFFT::help() {
    pprintf("-ifft performs an inverse dft on the current image, whose values are"
            " complex. The input and output are images with 2*c channels, where"
            " channel 2*i is the real part of the i\'th channel, and channel 2*i+1"
            " is the imaginary part of the i'th channel.\n"
            "\n"
            "Usage: ImageStack -load a.tga -fftcomplex -save freq.tga\n\n");
}

bool IFFT::test() {
    // tested by fft
    return true;
}

void IFFT::parse(vector<string> args) {
    assert(args.size() < 2, "-ifft takes zero or one argument\n");

    bool x = true, y = true, t = true;
    if (args.size() == 1) {
        x = y = t = false;
        for (size_t i = 0; i < args[0].size(); i++) {
            switch (args[0][i]) {
            case 'x':
                x = true;
                break;
            case 'y':
                y = true;
                break;
            case 't':
                t = true;
                break;
            default:
                panic("Unknown dimension: %c\n", args[0][i]);
                break;
            }
        }
    }

    apply(stack(0), x, y, t);
}


void IFFT::apply(Image im, bool x, bool y, bool t) {
    FFT::apply(im, x, y, t, true);
}

void FFTConvolve::help() {
    pprintf("-fftconvolve performs convolution in Fourier space. It is much faster"
            " than -convolve for large kernels. The two arguments are the boundary"
            " condition (zero, clamp, wrap, homogeneous) and the vector-vector"
            " multiplication used (inner, outer, elementwise). The defaults are wrap"
            " and outer respectively. See -convolve for a description of each"
            " option.\n"
            "\n"
            "Usage: ImageStack -load filter.tmp -load im.jpg -fftconvolve zero inner\n");
}

bool FFTConvolve::test() {
    Image im(17, 34, 8, 2);
    Image kernel(5, 7, 3, 2);
    Noise::apply(im, 0, 1);
    Noise::apply(kernel, 0, 1);
    // Make the kernel lopsided to check for flipping issues
    kernel(0, 0, 0, 0) = 17;

    Convolve::BoundaryCondition b[] = {Convolve::Clamp, Convolve::Wrap, Convolve::Zero, Convolve::Homogeneous};
    Multiply::Mode m[] = {Multiply::Outer, Multiply::Elementwise, Multiply::Inner};
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 3; j++) {
            Image fa = FFTConvolve::apply(im, kernel, b[i], m[j]);
            Image a = Convolve::apply(im, kernel, b[i], m[j]);
            if (!nearlyEqual(a, fa)) return false;
        }
    }

    return true;
}

void FFTConvolve::parse(vector<string> args) {
    Multiply::Mode m = Multiply::Outer;
    Convolve::BoundaryCondition b = Convolve::Wrap;

    if (args.size() > 0) {
        if (args[0] == "zero") { b = Convolve::Zero; }
        else if (args[0] == "homogeneous") { b = Convolve::Homogeneous; }
        else if (args[0] == "clamp") { b = Convolve::Clamp; }
        else if (args[0] == "wrap") { b = Convolve::Wrap; }
        else {
            panic("Unknown boundary condition: %s\n", args[0].c_str());
        }
    }

    if (args.size() > 1) {
        if (args[1] == "inner") { m = Multiply::Inner; }
        else if (args[1] == "outer") { m = Multiply::Outer; }
        else if (args[1] == "elementwise") { m = Multiply::Elementwise; }
        else {
            panic("Unknown vector-vector multiplication: %s\n", args[1].c_str());
        }
    }

    Image im = apply(stack(0), stack(1), b, m);
    pop();
    push(im);
}

Image FFTConvolve::apply(Image im, Image filter, Convolve::BoundaryCondition b, Multiply::Mode m) {

    Image out;
    if (m == Multiply::Inner) {
        assert(filter.channels % im.channels == 0 ||
               im.channels % filter.channels == 0,
               "To perform an inner or matrix product, the channel count "
               "of either the image or the filter must be a multiple of "
               "the channel count of the other.");
        if (im.channels < filter.channels) {
            out = Image(im.width, im.height, im.frames, filter.channels/im.channels);
            for (int i = 0; i < filter.channels; i++) {
                convolveSingle(im.channel(i % im.channels),
                               filter.channel(i),
                               out.channel(i / im.channels), b);
            }
        } else {
            out = Image(im.width, im.height, im.frames, im.channels/filter.channels);
            for (int i = 0; i < im.channels; i++) {
                convolveSingle(im.channel(i),
                               filter.channel(i % filter.channels),
                               out.channel(i / filter.channels), b);
            }
        }
    } else if (m == Multiply::Outer) {
        out = Image(im.width, im.height, im.frames, im.channels * filter.channels);
        for (int i = 0; i < im.channels; i++) {
            for (int j = 0; j < filter.channels; j++) {
                convolveSingle(im.channel(i),
                               filter.channel(j),
                               out.channel(i*filter.channels + j), b);
            }
        }
    } else if (m == Multiply::Elementwise) {
        assert(im.channels == filter.channels,
               "For element-wise multiplication, the image "
               "and filter must have the same number of channels.");
        out = Image(im.width, im.height, im.frames, im.channels);
        for (int i = 0; i < im.channels; i++) {
            convolveSingle(im.channel(i), filter.channel(i), out.channel(i), b);
        }

    } else {
        panic("Unknown multiplication mode");
    }
    return out;
}

void FFTConvolve::convolveSingle(Image im, Image filter, Image out, Convolve::BoundaryCondition b) {

    // Deal with the homogeneous case recursively. This is slightly
    // inefficient because we construct and transform the filter
    // twice, but it makes the code much simpler
    if (b == Convolve::Homogeneous) {
        Image result = apply(im, filter, Convolve::Zero, Multiply::Outer);
        Image weight(im.width, im.height, im.frames, 1);
        weight.set(1.0f);
        Image resultW = apply(weight, filter, Convolve::Zero, Multiply::Outer);
        out += Stats(filter).sum() * result / resultW;
        return;
    }

    assert(filter.width % 2 == 1 &&
           filter.height % 2 == 1 &&
           filter.frames % 2 == 1,
           "The filter must have odd dimensions\n");

    int xPad = filter.width/2;
    int yPad = filter.height/2;
    int tPad = filter.frames/2;

    if (b == Convolve::Wrap) {
        xPad = yPad = tPad = 0;
    }

    Image weightT;

    Image imT = Image(im.width+xPad*2, im.height+yPad*2, im.frames+tPad*2, 2);

    //printf("1\n"); fflush(stdout);
    // 1) Make the padded complex image
    if (b == Convolve::Clamp) {
        for (int t = 0; t < imT.frames; t++) {
            int st = clamp(t-tPad, 0, im.frames-1);
            for (int y = 0; y < imT.height; y++) {
                int sy = clamp(y-yPad, 0, im.height-1);
                for (int x = 0; x < imT.width; x++) {
                    int sx = clamp(x-xPad, 0, im.width-1);
                    imT(x, y, t, 0) = im(sx, sy, st, 0);
                }
            }
        }
    } else { // Zero or Wrap
        imT.region(xPad, yPad, tPad, 0,
                   im.width, im.height, im.frames, 1).set(im);
    }

    //printf("2\n"); fflush(stdout);
    // 2) Transform the padded image
    FFT::apply(imT);

    //printf("3\n"); fflush(stdout);
    // 3) Make a padded complex filter of the same size
    Image filterT(imT.width, imT.height, imT.frames, 2);
    for (int t = 0; t < filter.frames; t++) {
        int ft = t - filter.frames/2;
        if (ft < 0) ft += filterT.frames;
        for (int y = 0; y < filter.height; y++) {
            int fy = y - filter.height/2;
            if (fy < 0) fy += filterT.height;
            for (int x = 0; x < filter.width; x++) {
                int fx = x - filter.width/2;
                if (fx < 0) fx += filterT.width;
                filterT(fx, fy, ft, 0) = filter(x, y, t, 0);
            }
        }
    }

    //printf("4\n"); fflush(stdout);
    // 4) Transform the padded filter
    FFT::apply(filterT);

    //printf("5\n"); fflush(stdout);
    // 5) Multiply the two into a padded complex transformed result
    ComplexMultiply::apply(imT, filterT);

    //printf("6\n"); fflush(stdout);
    // 6) Inverse transorm the result
    IFFT::apply(imT);

    //printf("7\n"); fflush(stdout);
    // 7) Remove the padding, and convert back to real numbers
    out += imT.region(xPad, yPad, tPad, 0,
                      im.width, im.height, im.frames, 1);
}


void FFTPoisson::help() {
    printf("-fftpoisson computes an image from a gradient field in the same way as"
           " -poisson. It interprets the top image on the stack as the y gradient,"
           " and the next image as the x gradient. If a single argument is given,"
           " it uses that as a weight, and interprets the third image on the stack"
           " as a rough target output. The output of this operation will adhere to"
           " the target proportionally to the given weight.\n"
           "\n"
           "Usage: ImageStack -load gx.tmp -load gy.tmp -fftpoisson -display\n");
}

bool FFTPoisson::test() {
    Image a(233, 123, 1, 5);
    Noise::apply(a, 0, 1);
    Image dx = a.copy();
    Gradient::apply(dx, 'x');
    Image dy = a.copy();
    Gradient::apply(dy, 'y');
    Image b = FFTPoisson::apply(dx, dy, Image(), 0.00001);
    return nearlyEqual(a, b);
}

void FFTPoisson::parse(vector<string> args) {
    Image im;

    if (args.size() == 0) {
        im = apply(stack(1), stack(0), Image(), 0);
    } else if (args.size() == 1) {
        im = apply(stack(1), stack(0), stack(2), readFloat(args[0]));
    } else {
        panic("-fftpoisson takes zero or one arguments\n");
    }

    push(im);
}

// This implementation was based on code by Pravin Bhat and is
// available at:
// http://grail.cs.washington.edu/projects/screenedPoissonEq/ Bhat P.,
// Curless B., Cohen M., and Zitnick L. Fourier Analysis of the 2D
// Screened Poisson Equation for Gradient Domain Problems. European
// Conference on Computer Vision (ECCV) 2008.

// It was modified for ImageStack by Neeraj Agrawal and Ritvik Mudur,
// and further modified by Andrew Adams to change the boundary
// conditions expected on the gradient images (ImageStack uses zero
// boundary conditions on gradient images).

Image FFTPoisson::apply(Image dx, Image dy, Image target, float targetStrength) {

    assert(dx.width == dy.width &&
           dx.height == dy.height &&
           dx.frames == dy.frames &&
           dx.channels == dy.channels,
           "x gradient must be same size as y gradient\n");
    if (target.defined()) {
        assert(target.width == dx.width &&
               target.height == dx.height &&
               target.frames == dx.frames &&
               target.channels == dx.channels,
               "target image must have the same size as the gradient images\n");
    }

    Image fftBuff(dx.width, dx.height, 1, 1);

    //compute two 1D lookup tables for computing the DCT of a 2D Laplacian on the fly
    Image ftLapY(1, dx.height, 1, 1);
    Image ftLapX(dx.width, 1, 1, 1);

    for (int x = 0; x < dx.width; x++) {
        ftLapX(x, 0) = 2.0f * cos((M_PI * x) / (dx.width - 1));
    }

    for (int y = 0; y < dx.height; y++) {
        ftLapY(0, y) = -4.0f + (2.0f * cos((M_PI * y) / (dx.height - 1)));
    }

    // Create a DCT-I plan, which is its own inverse.
    fftwf_plan fftPlan;
    fftPlan = fftwf_plan_r2r_2d(dx.height, dx.width,
                                &fftBuff(0, 0), &fftBuff(0, 0),
                                FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE);

    Image out(dx.width, dx.height, dx.frames, dx.channels);

    for (int c = 0; c < dx.channels; c++) {
        for (int t = 0; t < dx.frames; t++) {

            float dcSum = 0.0f;

            // compute h_hat from u, gx, gy (see equation 48 in the paper), as well as the DC term of u's DCT.
            for (int y = 0; y < dx.height; y++) {
                for (int x = 0; x < dx.width; x++) {
                    // Compute DC term of u's DCT without computing the whole DCT.
                    float dcMult = 1.0f;
                    if ((x > 0) && (x < dx.width  - 1)) {
                        dcMult *= 2.0f;
                    }
                    if ((y > 0) && (y < dx.height - 1)) {
                        dcMult *= 2.0f;
                    }

                    if (target.defined()) {
                        dcSum += dcMult * target(x, y, t, c);
                    } else {
                        // try to read the dc term out of the double
                        // integral of the gradient fields
                        // instead. Works if the gradients were
                        // computed with a zero boundary condition.
                        dcSum += 2.0f*((dx.width-x)*dx(x, y, t, c) + (dy.height-y)*dy(x, y, t, c));
                    }


                    if (target.defined()) {
                        fftBuff(x, y) = targetStrength * target(x, y, t, c);
                    } else {
                        fftBuff(x, y) = 0;
                    }

                    // Subtract g^x_x and g^y_y, with boundary factor of -2.0 to account for boundary reflections implicit in the DCT
                    if (x == 0) {
                        fftBuff(x, y) -= (+2.0f * dx(x+1, y, t, c));
                    } else if (x == dx.width - 1) {
                        fftBuff(x, y) -= (-2.0f * dx(x, y, t, c));
                    } else {
                        fftBuff(x, y) -= (dx(x+1, y, t, c) - dx(x, y, t, c));
                    }

                    if (y == 0) {
                        fftBuff(x, y) -= (+2.0f * dy(x, y+1, t, c));
                    } else if (y == dx.height -1) {
                        fftBuff(x, y) -= (-2.0f * dy(x, y, t, c));
                    } else {
                        fftBuff(x, y) -= (dy(x, y+1, t, c) - dy(x, y, t, c));
                    }
                }
            }

            // transform h_hat to H_hat by taking the DCT of h_hat
            fftwf_execute(fftPlan);

            // compute F_hat using H_hat (see equation 29 in the paper)
            for (int y = 0; y < dx.height; y++) {
                for (int x = 0; x < dx.width; x++) {
                    float ftLapResponse = ftLapY(0, y) + ftLapX(x, 0);
                    fftBuff(x, y) /= (targetStrength - ftLapResponse);
                }
            }

            fftBuff(0, 0) = dcSum;

            // transform F_hat to f_hat by taking the inverse DCT of F_hat
            fftwf_execute(fftPlan);

            float fftMult = 1.0f / (4.0f * (dx.width-1) * (dx.height-1));

            for (int y = 0; y < dx.height; y++) {
                for (int x = 0; x < dx.width; x++) {
                    out(x, y, t, c) = fftBuff(x, y) * fftMult;
                }
            }

        }
    }

    fftwf_destroy_plan(fftPlan);

    return out;

}


#include "footer.h"
#endif

