#ifndef NO_FFTW
#include "main.h"
#include "Arithmetic.h"
#include "Complex.h"
#include "Deconvolution.h"
#include "DFT.h"
#include "File.h"
#include "GaussTransform.h"
#include "Geometry.h"
#include "KernelEstimation.h"
#include "Statistics.h"
#include "Filter.h"
#include "header.h"

#define FourierTransform(X) (FFT::apply(X, true, true, false))
#define InverseFourierTransform(X) (IFFT::apply(X, true, true, false))

void Deconvolve::help() {
    pprintf("-deconvolution will deconvolve an image with the kernel in the stack.\n"
            " This operation takes the name of the deconvolution method as a single\n"
            " argument, plus any optional arguments that the method may require.\n"
            " Currently supported are \"cho\" (Cho and Lee, 2009), \"shan\" \n"
            " (Shan et al, 2008), and \"levin\" - the simpler method with a"
            " Gaussian prior on gradients describe in (Levin 2007). The method"
            " \"levin\" takes an additional argument to specify the weight given to"
            " the prior.\n"
            "\n"
            "Usage: ImageStack -load blurred -load kernel -deconvolve cho\n");

}

bool Deconvolve::test() {
    // Make a lop-sided random kernel
    Image kernel(9, 9, 1, 1);
    Noise::apply(kernel, 0, 0.1);
    Noise::apply(kernel.region(0, 0, 0, 0,
                               5, 5, 1, 1), 0, 1);
    kernel /= Stats(kernel).sum();
    Image input = Downsample::apply(Load::apply("pics/dog1.jpg"), 2, 2, 1);
    Image blurry = Convolve::apply(input, kernel, Convolve::Zero);
    Noise::apply(blurry, -0.02, 0.02);
    Image shanResult = Deconvolve::applyShan2008(blurry, kernel);
    Image choResult = Deconvolve::applyCho2009(blurry, kernel);
    Image levinResult = Deconvolve::applyLevin2007(blurry, kernel, 0.02);
    /*
    Save::apply(kernel, "kernel.tmp");
    Save::apply(blurry, "blurry.tmp");
    Save::apply(shanResult, "shan.tmp");
    Save::apply(choResult, "cho.tmp");
    Save::apply(levinResult, "levin.tmp");
    */
    Stats shanStats(shanResult - input);
    Stats choStats(choResult - input);
    Stats levinStats(levinResult - input);
    printf("Shan:  %f %f\n"
           "Cho:   %f %f\n"
           "Levin: %f %f\n",
           shanStats.mean(), shanStats.variance(),
           choStats.mean(), choStats.variance(),
           levinStats.mean(), levinStats.variance());
    return (nearlyEqual(shanResult, input) &&
            nearlyEqual(choResult, input) &&
            nearlyEqual(levinResult, input));
}

void Deconvolve::parse(vector<string> args) {
    assert(args.size() >= 1, "-deconvolve takes at least one argument\n");
    Image kernel = stack(0);
    Image im = stack(1);
    if (args[0] == "cho") {
        assert(args.size() == 1, "-deconvolve cho takes no extra arguments\n");
        push(applyCho2009(im, kernel));
    } else if (args[0] == "shan") {
        assert(args.size() == 1, "-deconvolve shan takes no extra arguments\n");
        push(applyShan2008(im, kernel));
    } else if (args[0] == "levin") {
        assert(args.size() == 2, "-deconvolve levin takes one additional argument\n");
        push(applyLevin2007(im, kernel, readFloat(args[1])));
    } else {
        panic("Unknown method %s\n", args[0].c_str());
    }
}

Image Deconvolve::applyShan2008(Image B, Image K) {
    assert(K.channels == 1 && K.frames == 1,
           "The kernel must have one channel and one frame\n");
    assert(K.width % 2 == 1 && K.height % 2 == 1,
           "The kernel dimensions must be odd.\n");

    // Prepare constants and images.
    if (B.channels > 1 || B.frames > 1) {
        Image result(B.width, B.height, B.frames, B.channels);
        for (int c = 0; c < B.channels; c++) {
            for (int t = 0; t < B.frames; t++) {
                Image tmp = applyShan2008(B.channel(c).frame(t), K);
                result.channel(c).frame(t).set(tmp);
            }
        }
        return result;
    }
    Image B_large = applyPadding(B);
    Image K_large = KernelEstimation::enlargeKernel(K, B_large.width, B_large.height);
    Image smoothness_map;
    const int x_padding = (B_large.width - B.width) / 2;
    const int y_padding = (B_large.height - B.height) / 2;


    // Compute the smoothness map.
    {
        Image mean = B.copy();
        RectFilter::apply(mean, K.width | 1, K.height | 1, 1);
        Image variance = B * B;
        RectFilter::apply(variance, K.width | 1, K.height | 1, 1);
        smoothness_map = mean*mean - variance;
        Threshold::apply(smoothness_map, -25.0f / (256.f * 256.f));
        smoothness_map = Crop::apply(smoothness_map, -x_padding, -y_padding, 0, B_large.width, B_large.height, 1);
        Save::apply(smoothness_map, "smoothness_map.tmp");
    }

    // Prepare Fourier domain stuff.
    FourierTransform(K_large); // K_large = F(K).
    Image FK2 = K_large.copy(); // FK2 = F(K).
    ComplexConjugate::apply(K_large); // K_large = F(K)^T.
    ComplexMultiply::apply(FK2, K_large, false); // FK2 = |F(K)|^2.
    B_large = RealComplex::apply(B_large);
    FourierTransform(B_large);

    // sum w_i | K * (deriv_i L) - (deriv_i I) | ^ 2
    //   + gamma |Psi_x - deriv_x L|^2 + |Psi_y - deriv_y L|^2
    //        (Psi_x,Psi_y are redundant variables to follow deriv_x L, deriv_y L)
    //   + lambda_2 | Psi_x - deriv_x I|^2 + |Psi_y - deriv_y I|^2, masked by smoothness-Map
    //   + lambda_1 | non-linear prior on Psi_x, Psi_y |
    float lambda_1 = 0.1f, lambda_2 = 15.f;

    Image numerator_base(B_large.width, B_large.height, 1, 2);
    Image denominator_base(B_large.width, B_large.height, 1, 2);

    Image FDeriv[6];
    for (int i = 0; i <= 5; i++) {
        float w_i;
        FDeriv[i] = Image(B_large.width, B_large.height, 1, 2);
        switch (i) {
        case 0: // Original
            w_i = 50.f;
            FDeriv[i](0, 0) = 1.f; break;
        case 1: // dx
            w_i = 25.f;
            FDeriv[i](0, 0) = -1.f;
            FDeriv[i](1, 0) = 1.f; break;
        case 2: // dxx
            w_i = 12.5f;
            FDeriv[i](0, 0) = 1.f;
            FDeriv[i](1, 0) = -2.f;
            FDeriv[i](2, 0) = 1.f; break;
        case 3: // dy
            w_i = 25.f;
            FDeriv[i](0, 0) = -1.f;
            FDeriv[i](0, 1) = 1.f; break;
        case 4: // dyy
            w_i = 12.5f;
            FDeriv[i](0, 0) = 1.f;
            FDeriv[i](0, 1) = -2.f;
            FDeriv[i](0, 2) = 1.f; break;
        case 5: // dxy
            w_i = 12.5f;
            FDeriv[i](0, 0) = 1.f;
            FDeriv[i](1, 0) = -1;
            FDeriv[i](0, 1) = -1;
            FDeriv[i](1, 1) = 1; break;
        }
        FourierTransform(FDeriv[i]);
        Image tmp = FDeriv[i].copy();
        ComplexConjugate::apply(FDeriv[i]); // FDeriv[i] = F(deriv_i)^T
        ComplexMultiply::apply(tmp, FDeriv[i], false); // tmp = |F(deriv_i)|^2
        Image tmq = tmp.copy();
        ComplexMultiply::apply(tmp, FK2, false); // tmp = |F(K)|^2 |F(deriv_i)|^2
        ComplexMultiply::apply(tmq, K_large, false); // tmq = F(K)^T |F(deriv_i)|^2
        ComplexMultiply::apply(tmq, B_large, false); // tmq = F(K)^T |F(deriv_i)|^2 F(I)
        tmp *= w_i;
        tmq *= w_i;
        denominator_base += tmp;
        numerator_base += tmq;
    }

    Image dIdx = Convolve::apply(B_large, Crop::apply(FDeriv[1], -1, 0, 3, 1), Convolve::Wrap);
    Image dIdy = Convolve::apply(B_large, Crop::apply(FDeriv[3], 0, -1, 1, 3), Convolve::Wrap);
    Image L = B_large;
    Image FPsi_x, FPsi_y;
    Image Psi_x(B_large.width, B_large.height, 1, 2);
    Image Psi_y(B_large.width, B_large.height, 1, 2);
    float gamma = 2.0f;
    const int MAX_ITERATION = 1;

    // Non-linear prior for gradient (for 8-bit int pixels)
    float k = 2.7f;
    float a = 0.00061f;
    float b = 5.0f;
    float lt = 1.85263f;
    k *= 255.f; a *= 255.f * 255.f; lt /= 255.f; // adjustment for floating point

    for (int iterations = 1; iterations <= MAX_ITERATION; iterations++) {
        /******************************************/
        /* Optimize over Psi                      */
        /******************************************/
        // Fixing L makes the objective the following:
        //   gamma |Psi_x - deriv_x L|^2 + |Psi_y - deriv_y L|^2
        //   + lambda_2 | Psi_x - deriv_x I|^2 + |Psi_y - deriv_y I|^2, masked by smoothness-Map
        //   + lambda_1 | non-linear prior on Psi_x, Psi_y |
        // This can be done in the spatial domain entirely component-wise,
        // also independent for x and y.
        //  2 gamma (Psi_x - deriv_x L) + 2 lambda_2 (Psi_x - deriv_x I) .* mask
        //   + lambda_1 (non-linear prior on Psi_x)' = 0.
        float ans_x = 0, ans_y = 0, tmp, fscore = 0, tmpscore;
        bool fscore_valid;
        Image dLdx = Convolve::apply(L, Crop::apply(FDeriv[1], -1, 0, 3, 1), Convolve::Wrap);
        Image dLdy = Convolve::apply(L, Crop::apply(FDeriv[3], 0, -1, 1, 3), Convolve::Wrap);
        for (int y = 0; y < B_large.height; y++) {
            for (int x = 0; x < B_large.width; x++) {
                fscore_valid = false;
                // 1) Quadratic region:
                //  2 gamma (Psi_x - deriv_x L) + 2 lambda_2 (Psi_x - deriv_x I) .* mask
                //   + lambda_1 (-2a Psi_x) = 0.
                // Rearranging gives:
                //     (gamma + lambda_2 - a * lambda_1) psi_x = (gamma deriv_x L + lambda_2 deriv_x I .* mask)
                tmp = (gamma * dLdx(x, y) + lambda_2 * dIdx(x, y) * smoothness_map(x, y))
                      / (gamma + lambda_2 - a * lambda_1);
                tmpscore = gamma * (tmp - dLdx(x, y)) * (tmp - dLdx(x, y))
                           + lambda_2 * (tmp - dIdx(x, y)) * (tmp - dIdx(x, y)) * smoothness_map(x, y)
                           + lambda_1 * (-(a*tmp*tmp + b));
                if (fabs(tmp) > lt && (!fscore_valid || fscore > tmpscore)) {
                    fscore = tmpscore; ans_x = tmp; fscore_valid = true;
                }
                // 2) Positive linear region:
                //  2 gamma (Psi_x - deriv_x L) + 2 lambda_2 (Psi_x - deriv_x I) .* mask
                //   + lambda_1 (-k) = 0.
                // Rearranging gives:
                //     (gamma + lambda_2) psi_x = (gamma deriv_x L + lambda_2 deriv_x I .* mask + lambda_1 k)
                tmp = (gamma * dLdx(x, y) + lambda_2 * dIdx(x, y) * smoothness_map(x, y) + lambda_1 * k)
                      / (gamma + lambda_2);
                tmpscore = gamma * (tmp - dLdx(x, y)) * (tmp - dLdx(x, y))
                           + lambda_2 * (tmp - dIdx(x, y)) * (tmp - dIdx(x, y)) * smoothness_map(x, y)
                           + lambda_1 * (-k * tmp);
                if (tmp >= 0 && tmp <= lt && (!fscore_valid || fscore > tmpscore)) {
                    fscore = tmpscore; ans_x = tmp; fscore_valid = true;
                }
                // 3) Negative linear region:
                //  2 gamma (Psi_x - deriv_x L) + 2 lambda_2 (Psi_x - deriv_x I) .* mask
                //   + lambda_1 (k) = 0.
                // Rearranging gives:
                //     (gamma + lambda_2) psi_x = (gamma deriv_x L + lambda_2 deriv_x I .* mask - lambda_1 k)
                tmp = (gamma * dLdx(x, y) + lambda_2 * dIdx(x, y) * smoothness_map(x, y) - lambda_1 * k)
                      / (gamma + lambda_2);
                tmpscore = gamma * (tmp - dLdx(x, y)) * (tmp - dLdx(x, y))
                           + lambda_2 * (tmp - dIdx(x, y)) * (tmp - dIdx(x, y)) * smoothness_map(x, y)
                           + lambda_1 * (k * tmp);
                if (tmp >= 0 && tmp <= lt && (!fscore_valid || fscore > tmpscore)) {
                    fscore = tmpscore; ans_x = tmp; fscore_valid = true;
                }
                // 4) Zero
                tmp = 0.f;
                tmpscore = gamma * (tmp - dLdx(x, y)) * (tmp - dLdx(x, y))
                           + lambda_2 * (tmp - dIdx(x, y)) * (tmp - dIdx(x, y)) * smoothness_map(x, y);
                if (!fscore_valid || fscore > tmpscore) {
                    fscore = tmpscore; ans_x = tmp; fscore_valid = true;
                }
                // 5) lt
                tmp = lt;
                tmpscore = gamma * (tmp - dLdx(x, y)) * (tmp - dLdx(x, y))
                           + lambda_2 * (tmp - dIdx(x, y)) * (tmp - dIdx(x, y)) * smoothness_map(x, y)
                           - lambda_1 * (-k * lt);
                if (!fscore_valid || fscore > tmpscore) {
                    fscore = tmpscore; ans_x = tmp; fscore_valid = true;
                }
                // 6) -lt
                tmp = -lt;
                tmpscore = gamma * (tmp - dLdx(x, y)) * (tmp - dLdx(x, y))
                           + lambda_2 * (tmp - dIdx(x, y)) * (tmp - dIdx(x, y)) * smoothness_map(x, y)
                           - lambda_1 * (k * lt);
                if (!fscore_valid || fscore > tmpscore) {
                    fscore = tmpscore; ans_x = tmp; fscore_valid = true;
                }
                Psi_x(x, y) = ans_x;
            }
        }
        for (int y = 0; y < B_large.height; y++) {
            for (int x = 0; x < B_large.width; x++) {
                fscore_valid = false;
                // 1) Quadratic region:
                //  2 gamma (Psi_y - deriv_y L) + 2 lambda_2 (Psi_y - deriv_y I) .* mask
                //   + lambda_1 (-2a Psi_y) = 0.
                // Rearranging gives:
                //     (gamma + lambda_2 - a * lambda_1) psi_y = (gamma deriv_y L + lambda_2 deriv_y I .* mask)
                tmp = (gamma * dLdy(x, y) + lambda_2 * dIdy(x, y) * smoothness_map(x, y))
                      / (gamma + lambda_2 - a * lambda_1);
                tmpscore = gamma * (tmp - dLdy(x, y)) * (tmp - dLdy(x, y))
                           + lambda_2 * (tmp - dIdy(x, y)) * (tmp - dIdy(x, y)) * smoothness_map(x, y)
                           + lambda_1 * (-(a*tmp*tmp + b));
                if (fabs(tmp) > lt && (!fscore_valid || fscore > tmpscore)) {
                    fscore = tmpscore; ans_y = tmp; fscore_valid = true;
                }
                // 2) Positive linear region:
                //  2 gamma (Psi_y - deriv_y L) + 2 lambda_2 (Psi_y - deriv_y I) .* mask
                //   + lambda_1 (-k) = 0.
                // Rearranging gives:
                //     (gamma + lambda_2) psi_y = (gamma deriv_y L + lambda_2 deriv_y I .* mask + lambda_1 k)
                tmp = (gamma * dLdy(x, y) + lambda_2 * dIdy(x, y) * smoothness_map(x, y) + lambda_1 * k)
                      / (gamma + lambda_2);
                tmpscore = gamma * (tmp - dLdy(x, y)) * (tmp - dLdy(x, y))
                           + lambda_2 * (tmp - dIdy(x, y)) * (tmp - dIdy(x, y)) * smoothness_map(x, y)
                           + lambda_1 * (-k * tmp);
                if (tmp >= 0 && tmp <= lt && (!fscore_valid || fscore > tmpscore)) {
                    fscore = tmpscore; ans_y = tmp; fscore_valid = true;
                }
                // 3) Negative linear region:
                //  2 gamma (Psi_y - deriv_y L) + 2 lambda_2 (Psi_y - deriv_y I) .* mask
                //   + lambda_1 (k) = 0.
                // Rearranging gives:
                //     (gamma + lambda_2) psi_y = (gamma deriv_y L + lambda_2 deriv_y I .* mask - lambda_1 k)
                tmp = (gamma * dLdy(x, y) + lambda_2 * dIdy(x, y) * smoothness_map(x, y) - lambda_1 * k)
                      / (gamma + lambda_2);
                tmpscore = gamma * (tmp - dLdy(x, y)) * (tmp - dLdy(x, y))
                           + lambda_2 * (tmp - dIdy(x, y)) * (tmp - dIdy(x, y)) * smoothness_map(x, y)
                           + lambda_1 * (k * tmp);
                if (tmp >= 0 && tmp <= lt && (!fscore_valid || fscore > tmpscore)) {
                    fscore = tmpscore; ans_y = tmp; fscore_valid = true;
                }
                // 4) Zero
                tmp = 0.f;
                tmpscore = gamma * (tmp - dLdy(x, y)) * (tmp - dLdy(x, y))
                           + lambda_2 * (tmp - dIdy(x, y)) * (tmp - dIdy(x, y)) * smoothness_map(x, y);
                if (!fscore_valid || fscore > tmpscore) {
                    fscore = tmpscore; ans_y = tmp; fscore_valid = true;
                }
                // 5) lt
                tmp = lt;
                tmpscore = gamma * (tmp - dLdy(x, y)) * (tmp - dLdy(x, y))
                           + lambda_2 * (tmp - dIdy(x, y)) * (tmp - dIdy(x, y)) * smoothness_map(x, y)
                           - lambda_1 * (-k * lt);
                if (!fscore_valid || fscore > tmpscore) {
                    fscore = tmpscore; ans_y = tmp; fscore_valid = true;
                }
                // 6) -lt
                tmp = -lt;
                tmpscore = gamma * (tmp - dLdy(x, y)) * (tmp - dLdy(x, y))
                           + lambda_2 * (tmp - dIdy(x, y)) * (tmp - dIdy(x, y)) * smoothness_map(x, y)
                           - lambda_1 * (k * lt);
                if (!fscore_valid || fscore > tmpscore) {
                    fscore = tmpscore; ans_y = tmp; fscore_valid = true;
                }
                Psi_y(x, y) = ans_y;
            }
        }

        FPsi_x = Psi_x.copy(); FourierTransform(FPsi_x);
        FPsi_y = Psi_y.copy(); FourierTransform(FPsi_y);
        /******************************************/
        /* Optimize over L                        */
        /******************************************/
        // Fix Psi. Then the gradient of the objective becomes, in the Fourier domain:
        //    sum w_i F(K)^T F(deriv_i)^T (F(K) F(deriv_i) F(L) - F(deriv_i) F(I))
        //           + gamma F(deriv_x)^T (F(deriv_x) F(L) - F(Psi_x)) + ... = 0
        // Solving this yields F(L) = N/D,  where
        //   N = sum w_i F(K)^T |F(deriv_i)|^2 F(I) + gamma (F(deriv_x)^T F(Psi_x) +  ... )
        //   D = sum w_i |F(K)|^2 |F(deriv_i)|^2  + gamma |F(deriv_x)|^2 + |F(deriv_y)|^2
        // Note that the first term of N and D are independent of L or Psi or gamma.
        Image denominator = denominator_base.copy();
        Image numerator = numerator_base.copy();
        // Append the variable terms.
        for (int y = 0; y < B_large.height; y++)
            for (int x = 0; x < B_large.width; x++) {
                denominator(x, y) += gamma * (FDeriv[1](x, y, 0) * FDeriv[1](x, y, 0) +
                                              FDeriv[1](x, y, 1) * FDeriv[1](x, y, 1));
                denominator(x, y) += gamma * (FDeriv[3](x, y, 0) * FDeriv[3](x, y, 0) +
                                              FDeriv[3](x, y, 1) * FDeriv[3](x, y, 1));
                numerator(x, y) += gamma * (FDeriv[1](x, y, 0) * FPsi_x(x, y, 0) +
                                            FDeriv[1](x, y, 1) * FPsi_x(x, y, 1));
                numerator(x, y) += gamma * (FDeriv[3](x, y, 0) * FPsi_y(x, y, 0) +
                                            FDeriv[3](x, y, 1) * FPsi_y(x, y, 1));
            }
        ComplexDivide::apply(numerator, denominator, false);
        InverseFourierTransform(numerator);
        L = ComplexReal::apply(numerator);

        /*
        char filename_c[20];
        sprintf(filename_c, "output%02d.tmp", iterations);
        std::string filename(filename_c);
        FileTMP::save(L, filename, "float");
        */

        /******************************************/
        /* Bookkeeping                            */
        /******************************************/
        lambda_1 /= 1.2f; lambda_2 /= 1.5f; gamma *= 2.f;
    }
    // Crop L as before.
    L = Crop::apply(L, x_padding, y_padding, B.width, B.height);
    return L;
}

/*
 * A poor man's version of "Reducing Boundary Artifacts in Image Deconvolution" (ICCP 2008)
 * by Liu and Jia.
 */
Image Deconvolve::applyPadding(Image B) {
    // Calculate the margin size.
    int alpha = 1;
    if (B.width / 3 < alpha) alpha = B.width / 3;
    if (B.height / 3 < alpha) alpha = B.height / 3;
    int x_padding = B.width / 2;
    int y_padding = B.height / 2;
    if (x_padding < alpha * 3) x_padding = alpha * 3;
    if (y_padding < alpha * 3) y_padding = alpha * 3;

    // Prepare the enlarged canvas.
    vector<float> prev(B.channels);
    Image ret = Crop::apply(B, -x_padding, -y_padding, 0,
                            B.width+x_padding*2, B.height+y_padding*2, B.frames);
    for (int t = 0; t < B.frames; t++) {
        // Populate the top 'A' region.
        for (int y = 0; y < alpha; y++) {
            for (int c = 0; c < B.channels; c++) {
                for (int x = 0; x < B.width; x++) {
                    ret(x+x_padding, y, t, c) = ret(x+x_padding, y - alpha + B.height + y_padding, t, c);
                    ret(x+x_padding, y_padding - alpha + y, t, c) = ret(x+x_padding, y + y_padding, t, c);
                }
            }
        }
        for (int y = alpha; y < y_padding - alpha; y++) {
            // interpolate towards the bottom boundary.
            float weight = 1.f / (y_padding - alpha - (y-1));
            for (int x = x_padding; x < x_padding + B.width; x++)
                for (int c = 0; c < B.channels; c++)
                    ret(x, y, t, c) = ret(x, y-1, t, c) * (1.f - weight) + ret(x, y_padding - alpha, t, c) * weight;
            // Blur with neighbors, more increasingly at the center.
            for (int c = 0; c < B.channels; c++)
                prev[c] = ret(x_padding, y, t, c);
            float wing = 0.1f + 0.2f * (1.f - fabs(y_padding * 0.5f - y) / (y_padding * 0.5f));
            float center = 1.f - wing * 2.f;
            for (int x = x_padding; x < x_padding + B.width - 1; x++) {
                for (int c = 0; c < B.channels; c++) {
                    float tmp = ret(x, y, t, c);
                    ret(x, y, t, c) = prev[c] * wing + ret(x+1, y, t, c) * wing + tmp * center;
                    prev[c] = tmp;
                }
            }
        }
        // Populate the bottom 'A' region
        for (int y = 0; y < y_padding; y++) {
            for (int c = 0; c < B.channels; c++) {
                for (int x = 0; x < B.width; x++) {
                    ret(x+x_padding, y+B.height+y_padding, t, c) = ret(x+x_padding, y, t, c);
                }
            }
        }

        // Populate the left 'C-B-C' region
        for (int y = 0; y < B.height + y_padding * 2; y++) {
            for (int x = 0; x < alpha; x++) {
                for (int c = 0; c < B.channels; c++) {
                    ret(x, y, t, c) = ret(B.width + x_padding - alpha + x, y, t, c);
                    ret(x_padding - alpha + x, y, t, c) = ret(x_padding + x, y, t, c);
                }
            }
        }

        for (int x = alpha; x < x_padding - alpha; x++) {
            // interpolate towards the right boundary.
            float weight = 1.f / (x_padding - alpha - (x-1));
            for (int y = 0; y < B.height + y_padding * 2; y++) {
                for (int c = 0; c < B.channels; c++) {
                    ret(x, y, t, c) = ret(x-1, y, t, c) * (1.f - weight) + ret(x_padding - alpha, y, t, c) * weight;
                }
            }
            // Blur with neighbors, more increasingly at the center.
            for (int c = 0; c < B.channels; c++) {
                prev[c] = ret(x, 0, t, c);
            }
            float wing = 0.1f + 0.2f * (1.f - fabs(x_padding * 0.5f - x) / (x_padding * 0.5f));
            float center = 1.f - wing * 2.f;
            for (int y = 0; y < B.height + y_padding * 2 - 1; y++) {
                for (int c = 0; c < B.channels; c++) {
                    float tmp = ret(x, y, t, c);
                    ret(x, y, t, c) = prev[c] * wing + ret(x, y+1, t, c) * wing + tmp * center;
                    prev[c] = tmp;
                }
            }
        }
        // Populate the right 'C-B-C' region
        for (int y = 0; y < B.height + y_padding * 2; y++) {
            for (int c = 0; c < B.channels; c++) {
                for (int x = 0; x < x_padding; x++) {
                    ret(x + B.width + x_padding, y, t, c) = ret(x, y, t, c);
                }
            }
        }
    }
    ret = Crop::apply(ret, x_padding/2, y_padding/2, 0, B.width + x_padding, B.height + y_padding, B.frames);
    return ret;
}

Image Deconvolve::applyCho2009(Image blurred, Image kernel) {
    assert(kernel.width % 2 == 1 && kernel.height % 2 ==1,
           "The kernel dimensions must be odd.\n");
    assert(kernel.channels == 1 && kernel.frames == 1 && blurred.frames == 1,
           "The kernel must be single-channel, and both the kernel and blurred\n"
           "image must be single-framed.\n");

    // omega_* = 50 / (2^q) where q is the order of the derivative, alpha = 0.1
    // want to minimize w.r.t. L:

    // sum w_i | K * (deriv_i L) - (deriv B) | ^ 2 + alpha | grad L | ^ 2
    // In Fourier domain,
    //  sum w_i |F(K) F(deriv_i L) - F(deriv_i B)|^2 + alpha |F(grad L)|^2
    // sum w_i |F(K).* F(deriv_i).*F(L) - F(deriv_i).*F(B)|^2 + alpha ( |F(dx).*F(L)|^2 + |F(dy).*F(L)|^2)
    // Differentiating w.r.t. each element of F(L) and setting to zero, we get
    // F(L) = F(K)^T F(B) sum_i w_i |F(deriv_i)|^2  divided by
    //          |F(K)|^2 sum_i w_i |F(deriv_i)|^2  + alpha (|F(dx)|^2+|F(dy)|^2)

    Image B  = applyPadding(blurred);
    //FileTMP::save(B, std::string("padded.tmp"), "float");

    float alpha = 1.f; // TODO
    Image FK = KernelEstimation::enlargeKernel(kernel, B.width, B.height);
    Image FB = RealComplex::apply(Transpose::apply(B, 'c', 't'));
    FFT::apply(FK, true, true, false);
    FFT::apply(FB, true, true, false);
    Image FK2 = FK.copy();
    ComplexMultiply::apply(FK2, FK, true);
    Image SumDeriv(B.width, B.height, 1, 2);
    Image SumGrad(B.width, B.height, 1, 2);
    for (int i = 0; i <= 5; i++) {
        float w_i;
        Image FDeriv(B.width, B.height, 1, 2);
        switch (i) {
        case 0: // Original
            w_i = 50.f;
            FDeriv(0, 0) = 1.f; break;
        case 1: // dx
            w_i = 25.f;
            FDeriv(0, 0) = -1.f;
            FDeriv(1, 0) = 1.f; break;
        case 2: // dxx
            w_i = 12.5f;
            FDeriv(0, 0) = 1.f;
            FDeriv(1, 0) = -2.f;
            FDeriv(2, 0) = 1.f; break;
        case 3: // dy
            w_i = 25.f;
            FDeriv(0, 0) = -1.f;
            FDeriv(0, 1) = 1.f; break;
        case 4: // dyy
            w_i = 12.5f;
            FDeriv(0, 0) = 1.f;
            FDeriv(0, 1) = -2.f;
            FDeriv(0, 2) = 1.f; break;
        case 5: // dxy
            w_i = 12.5f;
            FDeriv(0, 0) = 1.f;
            FDeriv(1, 0) = -1;
            FDeriv(0, 1) = -1;
            FDeriv(1, 1) = 1; break;
        }
        FFT::apply(FDeriv, true, true, false);
        Image FDeriv2 = FDeriv.copy();
        ComplexMultiply::apply(FDeriv2, FDeriv, true);
        if (i == 1 || i == 3) {
            SumGrad += FDeriv2;
        }
        FDeriv2 *= w_i;
        SumDeriv += FDeriv2;
    }
    SumGrad *= alpha;
    // Recall the following:
    // F(L) = F(K)^T F(B) sum_i w_i |F(deriv_i)|^2  divided by
    //          |F(K)|^2 sum_i w_i |F(deriv_i)|^2  + alpha (|F(dx)|^2+|F(dy)|^2)
    // In our diction, we have
    // FK^T FB SumDeriv / (FK2 SumDeriv + SumGrad)
    ComplexConjugate::apply(FK);
    ComplexMultiply::apply(FK, SumDeriv, false);
    ComplexMultiply::apply(FK2, SumDeriv, false);
    FK2 += SumGrad;
    ComplexDivide::apply(FK, FK2, false); // FK contains the running result.
    for (int t = 0; t < B.channels; t++) {
        for (int y = 0; y < B.height; y++) {
            for (int x = 0; x < B.width; x++) {
                float tmp = FB(x, y, t, 0) * FK(x, y, 0) - FB(x, y, t, 1) * FK(x, y, 1);
                FB(x, y, t, 1) = FB(x, y, t, 0) * FK(x, y, 1) + FB(x, y, t, 1) * FK(x, y, 0);
                FB(x, y, t, 0) = tmp;
            }
        }
    }
    IFFT::apply(FB, true, true, false);
    const int x_padding = (B.width - blurred.width) / 2;
    const int y_padding = (B.height - blurred.height) / 2;
    return Crop::apply(Transpose::apply(ComplexReal::apply(FB), 'c', 't'),
                       x_padding, y_padding, 0, blurred.width, blurred.height, blurred.frames);
}

Image Deconvolve::applyLevin2007(Image blurred, Image kernel, float weight) {
    assert(kernel.width % 2 == 1 && kernel.height % 2 ==1,
           "The kernel dimensions must be odd.\n");
    assert(kernel.channels == 1 && kernel.frames == 1 && blurred.frames == 1,
           "The kernel must be single-channel, and both the kernel and blurred\n"
           "image must be single-framed.\n");

    Image padded = applyPadding(blurred);

    // sum of second derivatives filter
    Image fft_g(padded.width, padded.height, 1, 2);
    fft_g(0, 0) = weight;
    fft_g(padded.width-1, 0) = -weight*0.25;
    fft_g(0, padded.height-1) = -weight*0.25;
    fft_g(1, 0) = -weight*0.25;
    fft_g(0, 1) = -weight*0.25;
    FFT::apply(fft_g);

    Image fft_im = RealComplex::apply(padded);
    FFT::apply(fft_im);

    Image fft_kernel(padded.width, padded.height, 1, 2);
    for (int y = 0; y < kernel.height; y++) {
        int fy = y - kernel.height/2;
        if (fy < 0) { fy += fft_kernel.height; }
        for (int x = 0; x < kernel.width; x++) {
            for (int c = 0; c < kernel.channels; c++) {
                int fx = x - kernel.width/2;
                if (fx < 0) { fx += fft_kernel.width; }
                fft_kernel(fx, fy, 2*c) = kernel(x, y, c);
            }
        }
    }
    FFT::apply(fft_kernel);

    ComplexMultiply::apply(fft_im, fft_kernel, true);
    ComplexMultiply::apply(fft_kernel, fft_kernel, true);

    fft_kernel += fft_g;
    ComplexDivide::apply(fft_im, fft_kernel, false);


    const int x_pad = (padded.width - blurred.width)/2;
    const int y_pad = (padded.height - blurred.height)/2;

    IFFT::apply(fft_im);
    fft_im = fft_im.region(x_pad, y_pad, 0, 0,
                           blurred.width, blurred.height,
                           blurred.frames, fft_im.channels);
    return ComplexReal::apply(fft_im);
}

#include "footer.h"

#endif
