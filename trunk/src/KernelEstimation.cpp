#ifndef NO_FFTW
#include "main.h"
#include "Arithmetic.h"
#include "Color.h"
#include "Complex.h"
#include "DFT.h"
#include "Deconvolution.h"
#include "File.h"
#include "Geometry.h"
#include "KernelEstimation.h"
#include "header.h"

#define FourierTransform(X) (FFT::apply(X, true, true, false))
#define InverseFourierTransform(X) (IFFT::apply(X, true, true, false))

void KernelEstimation::help() {
    pprintf("-kernelestimation will compute an estimated blur kernel from a blurry"
            " image using the algorithm described in Cho and Lee, 2009. It takes an"
            " optional argument that specifies the kernel size. The default is 25.\n"
            "Usage: ImageStack -load blurred -kernelestimation 25 -deconvolve cho\n");
}

void KernelEstimation::parse(vector<string> args) {
    assert(args.size() <= 1, "-kernelestimation takes at most one argument.\n");
    int kernel_size = args.size() == 1 ? readInt(args[0]) : 25;
    Image im = apply(stack(0), kernel_size);
    push(im);
}

/*
 * Normalizes the window to have a sum of 1, provided that
 * the window has one channel and one frame, and a nonzero sum.
 */
void KernelEstimation::NormalizeSum(Window im) {
    assert(im.channels == 1 && im.frames == 1,
           "The image to be normalized must have one channel and one frame!\n");
    double sum = 0.0;
    for (int y = 0; y < im.height; y++)
        for (int x = 0; x < im.width; x++)
            sum += im(x, y)[0];
    Scale::apply(im, 1.0 / sum);
}

Image KernelEstimation::EnlargeKernel(Window im, int w, int h) {
    Image ret(w, h, 1, 2);
    for (int y = 0; y < im.height; y++) {
        int new_y = (y - (im.height >> 1) + h) % h;
        int new_x =  - (im.width >> 1) + w;
        for (int x = 0; x < im.width; x++, new_x++) {
            if (new_x == w) new_x = 0;
            ret(new_x, new_y)[0] = im(x, y)[0];
        }
    }
    return ret;
}

Image KernelEstimation::ContractKernel(Window im, int size) {
    Image ret(size, size, 1, 1);
    for (int y = 0; y < size; y++) {
        int y_old = (y - (size >> 1) + im.height) % im.height;
        int x_old = ( - (size >> 1) + im.width) % im.width;
        for (int x = 0; x < size; x++, x_old++) {
            if (x_old >= im.width) x_old = 0;
            ret(x, y)[0] = im(x_old, y_old)[0];
        }
    }
    return ret;
}

/*
 * Applies the shock filter of Osher and Rudin 1990, as described by Cho and Lee 2009.
 */
void KernelEstimation::ShockFilterIteration(Window im, float dt) {
    Image input(im.width+2, im.height+2, im.frames, im.channels);
    Paste::apply(input, im, 0, 0, 0, 0, im.width, im.height);
    for (int t = 0; t < im.frames; t++) {
        for (int y = 1; y < im.height - 1; y++) {
            for (int x = 1; x < im.width - 1; x++) {
                for (int c = 0; c < im.channels; c++) {
                    float dx = input(x, y, t)[c] - input(x+1, y, t)[c];
                    float dy = input(x, y, t)[c] - input(x, y+1, t)[c];
                    float laplacian = input(x, y-1, t)[c] + input(x, y+1, t)[c] +
                        input(x-1, y, t)[c] + input(x+1, y, t)[c] - 4 * input(x, y, t)[c];
                    im(x, y, t)[c] -= sqrtf(dx*dx+dy*dy) * (laplacian > 0.f ? 1.f : -1.f) * dt;
                }
            }
        }
    }
}

/*
 * Applies the bilateral filter with the parameters described in Cho and Lee 2009.
 * The support is a 5x5 window, with spatial standard deviation at 2, and the range
 * standard deviation specified.
 */
void KernelEstimation::BilateralFilterIteration(Window im, float sigma_r) {
    float sigma_s = 2.0f;
    float gaussian[5][5];
    for (int y = 0; y < 5; y++)
        for (int x = 0; x < 5; x++)
            gaussian[y][x] = expf(-(float)((x-2)*(x-2)+(y-2)*(y-2)) / (2.0f * sigma_s * sigma_s));
    Image input(im);
    for (int t = 0; t < im.frames; t++)
        for (int y = 0; y < im.height; y++)
            for (int x = 0; x < im.width; x++)
                for (int c = 0; c < im.channels; c++) {
                    float total_weight = 0.f, pixel = 0.f;
                    for (int y2 = (y>2?y-2:0); y2 <= y + 2 && y2 < im.height; y2++) {
                        for (int x2 = (x>2?x-2:0); x2 <= x + 2 && x2 < im.width; x2++) {
                            float color_dist = input(x2,y2,t)[c] - input(x,y,t)[c];
                            float weight = expf(-color_dist * color_dist / (2.0f * sigma_r * sigma_r)) * gaussian[y2-y+2][x2-x+2];
                            total_weight += weight;
                            pixel += weight * input(x2,y2,t)[c];
                        }
                    }
                    im(x,y,t)[c] = pixel / total_weight;
                }
}

/*
 * Computes the component-wise multiplication between the two images and reports the sum.
 */
float KernelEstimation::DotProduct(Window im1, Window im2, int channel, int frame) {
    assert(im1.width == im2.width && im1.height == im2.height &&
           im1.channels > channel && im2.channels > channel &&
           im1.frames > frame && im2.frames > frame,
           "Either the two images have mismatched dimensions, or not enough "
           "channels or frames.\n");
    double ret = 0.f;
    for (int y = 0; y < im1.height; y++)
        for (int x = 0; x < im1.width; x++)
            ret += im1(x, y, frame)[channel] * im2(x, y, frame)[channel];
    return (float)ret;
}

Image KernelEstimation::BilinearResample(Image im, int w, int h) {
    Image ret(w, h, im.frames, im.channels);
    float x_factor = ((float)im.width-1) / (w-1);
    for (int t = 0; t < im.frames; t++) {
        for (int y = 0; y < h - 1; y++) {
            float y_old = ((float)y) / (h-1) * (im.height-1);
            int y_old_int = (int)y_old;
            float y_old_frac = y_old - y_old_int;
            for (int x = 0; x < w - 1; x++)
                im.sample2DLinear(x * x_factor, y_old, t, ret(x, y, t));
            for (int c = 0; c < im.channels; c++) {
                ret(w-1, y, t)[c] = im(im.width-1, y_old_int, t)[0] * (1 - y_old_frac)
                    + im(im.width-1, y_old_int + 1, t)[0] * (y_old_frac);
            }
        }
        for (int x = 0; x < w - 1; x++) {
            float x_old = x * x_factor;
            int x_old_int = (int)x_old;
            float x_old_frac = x_old - x_old_int;
            for (int c = 0; c < im.channels; c++) {
                ret(x, h-1, t)[c] = im(x_old_int, im.height-1, t)[0] * (1 - x_old_frac)
                    + im(x_old_int+1, im.height-1, t)[0] * (x_old_frac);
            }
        }
        for (int c = 0; c < im.channels; c++)
            ret(w-1,h-1,t)[c] = im(im.width-1,im.height-1,t)[c];
    }
    return ret;
}

Image KernelEstimation::apply(Window B, int kernel_size) {

    /******************************* Parameter check */

    assert((kernel_size % 2) == 1, "The desired kernel size must be odd.\n");

    assert((B.channels == 1 || B.channels == 3) && B.frames == 1,
           "The blurred image must have 1 or 3 channels, and 1 frame.\n");

    /******************************* Initialize Constants */
    float dt = 1.f;
    float sigma_r = 0.5f;
    float gradient_threshold = 0.f;
    const int CG_ITERATIONS = 10;
    int primes[] = {2, 3, 5, 7};

    /******************************* Compute Kernel Sizes */
    std::vector<int> kernel_scale;
    int tmp = kernel_size;
    kernel_scale.push_back(tmp);
    while (tmp != 3) {
        tmp = ((int)((float)tmp / 1.6f) / 2) * 2 + 1;
        kernel_scale.push_back(tmp);
    }

    /******************************* Declare Local Variables */
    Image Bgray = (B.channels == 3) ? ColorConvert::apply(B, "rgb", "y") : Image(B);
    Image Blurry;
    Image guess;
    Image K(3, 3, 1, 1);
    K(1,1)[0] = 0.4;
    K(0,1)[0] = K(2,1)[0] = K(1,0)[0] = K(1,2)[0] = 0.15;

    float tic = currentTime(), toc;
    /******************************* Coarse-to-Fine Loop */
    for (unsigned int iteration = 1; iteration <= kernel_scale.size(); iteration++) {
        // Setup.
        int m = kernel_scale[kernel_scale.size()-iteration];
        int newwidth = ((float)m) / kernel_size * B.width;
        int newheight = ((float)m) / kernel_size * B.height;
        int padded_width, padded_height;
        for (int i = 0; i < 4; i++) {
            int s;
            for (s = 1; s < newwidth + m - 1; s *= primes[i]);
            if (s < padded_width || i == 0) padded_width = s;
            for (s = 1; s < newheight + m - 1; s *= primes[i]);
            if (s < padded_height || i == 0) padded_height = s;
        }
        printf("Iteration %d/%d (%dx%d, %dx%d)\n", iteration,
               (int)kernel_scale.size(), newwidth, newheight, m, m);

        // Generate the important images at the current scale.
        Blurry = BilinearResample(Bgray, newwidth, newheight);

        //    if (iteration <= 3)
        //      guess = ColorConvert::apply(Load::apply("DeblurTest/input3.tmp"), "rgb", "y");

        guess = BilinearResample(iteration == 1 ? Bgray : guess, newwidth, newheight);

        K = BilinearResample(K, m, m);
        NormalizeSum(K);
        toc = currentTime(); printf(" Setup     : %.3f sec\n", toc - tic); tic = toc;

        /**************************************************************/
        /* PREDICTION                                                 */
        /**************************************************************/
        BilateralFilterIteration(guess, sigma_r);
        ShockFilterIteration(guess, dt);

        /*char filename_c[20];
        sprintf(filename_c, "prediction%d.tmp", iteration);
        FileTMP::save(guess, std::string(filename_c), "float");*/

        // Compute gradient and generate an angular histogram.
        if (iteration == 1) {
            // Compute the gradient histogram.
            unsigned int GRAD_CAP = min((newwidth-1) * (newheight-1), 100000);
            float probability_cutoff = GRAD_CAP / (newwidth-1) * (newheight-1);
            std::vector<float> GRAD_HIST[4];
            for (int i = 0; i < 4; i++) {
                GRAD_HIST[i].clear();
                GRAD_HIST[i].reserve(GRAD_CAP);
            }
            for (int y = 0; y < newheight - 1; y++)
                for (int x = 0; x < newwidth - 1; x++) {
                    float dx = guess(x+1, y)[0] - guess(x, y)[0];
                    float dy = guess(x, y+1)[0] - guess(x, y)[0];
                    float angle = atan2(dy, dx);
                    int bin = (int)(floor(angle / (M_PI * 0.25f)) + 4) % 4;
                    if (GRAD_HIST[bin].size() < GRAD_CAP && rand() / RAND_MAX < probability_cutoff) {
                        GRAD_HIST[bin].push_back(dx*dx+dy*dy);
                    }
                }
            // Find the cutoff.
            float cutoff = 0.f;
            float inclusion_count = m * 4.0f; // m * r in the paper, and r=2.
            for (int i = 0; i < 4; i++) {
                std::sort(GRAD_HIST[i].begin(), GRAD_HIST[i].end());
                int cutoff_index = 0;
                cutoff_index = GRAD_HIST[i].size() - 1 - inclusion_count;
                if (cutoff_index >= 0 && cutoff < GRAD_HIST[i][cutoff_index])
                    cutoff = GRAD_HIST[i][cutoff_index];
            }
            gradient_threshold = cutoff * 0.1f;
        } else {
            gradient_threshold *= 0.9f * 0.9f;
        }
        // Generate the gradient channels.
        Image Px(padded_width, padded_height, 1, 2);
        Image Py(padded_width, padded_height, 1, 2);
        Image Bx(padded_width, padded_height, 1, 2);
        Image By(padded_width, padded_height, 1, 2);
        int xoffset = (padded_width - newwidth) >> 1;
        int yoffset = (padded_height - newheight) >> 1;
        for (int y = 0; y < newheight - 1; y++) {
            for (int x = 0; x < newwidth - 1; x++) {
                float dx = guess(x+1, y)[0] - guess(x, y)[0];
                float dy = guess(x, y+1)[0] - guess(x, y)[0];
                if (dx*dx+dy*dy > gradient_threshold) {
                    Px(x + xoffset, y + yoffset)[0] = dx; Py(x + xoffset, y + yoffset)[0] = dy;
                }
                Bx(x + xoffset, y + yoffset)[0] = Blurry(x+1, y)[0] - Blurry(x, y)[0];
                By(x + xoffset, y + yoffset)[0] = Blurry(x, y+1)[0] - Blurry(x, y)[0];
            }
        }
        // Update constants for next iteration.
        dt *= 0.9f;
        sigma_r *= 0.9f;    
        toc = currentTime(); printf(" Prediction: %.3f sec\n", toc - tic); tic = toc;
        /*
        sprintf(filename_c, "edges%d.tmp", iteration);
        FileTMP::save(Adjoin::apply(Px, Py, 't'), std::string(filename_c), "float");
        */

        /**************************************************************/
        /* KERNEL ESTIMATION                                          */
        /**************************************************************/
        // Build the gradient images.
        float beta = 1.f;
        Image dxPx(padded_width, padded_height, 1, 2);    
        Image dyPy(padded_width, padded_height, 1, 2);    
        Image dxyPxy(padded_width, padded_height, 1, 2);    
        Image dxBx(padded_width, padded_height, 1, 2);    
        Image dyBy(padded_width, padded_height, 1, 2);    
        Image dxyBxy(padded_width, padded_height, 1, 2);    
        for (int y = yoffset; y < newheight - 2 + yoffset; y++) {
            for (int x = xoffset; x < newwidth - 2 + xoffset; x++) {
                dxPx(x, y)[0] = Px(x+1, y)[0] - Px(x, y)[0];
                dyPy(x, y)[0] = Py(x, y+1)[0] - Py(x, y)[0];
                dxyPxy(x, y)[0] = (Px(x, y+1)[0] - Px(x, y)[0] + Py(x+1, y)[0] - Py(x, y)[0]) * 0.5f;
                dxBx(x, y)[0] = Bx(x+1, y)[0] - Bx(x, y)[0];
                dyBy(x, y)[0] = By(x, y+1)[0] - By(x, y)[0];
                dxyBxy(x, y)[0] = (Bx(x, y+1)[0] - Bx(x, y)[0] + By(x+1, y)[0] - By(x, y)[0]) * 0.5f;
            }
        }
        FourierTransform(Px); FourierTransform(Py);
        FourierTransform(dxPx); FourierTransform(dyPy); FourierTransform(dxyPxy);
        FourierTransform(Bx); FourierTransform(By);
        FourierTransform(dxBx); FourierTransform(dyBy); FourierTransform(dxyBxy);
        // Need to write conjugate gradient minimizing the following objective:
        // f(K) = sum_i w_i|A_iK-B_i|^2 + beta |K|^2
        //   where B_i is the i-th derivative of blurry image
        //          and A_i is the i-th derivative of the latent image
        //           expressed in a matrix form (as coefficients to K).
        //
        // In quadratic form, this comes out to be:
        // f(K) = 1/2 K^T 2(sum_i w_i A_i^TA_i + beta 1) K  - 2 (sum w_i A_i^TB_i) K
        //      = 1/2 K^T CoeffA K - CoeffB K
        // Initialize a guess for the kernel.

        // This requires us to precompute, among other things,
        // "CoeffA" = 2.0 * sum_i w_i F{A_i^T} .* F{A_i}) + F{beta} /// in fourier domain!
        // "CoeffB" = 2.0 * sum_i w_i F^{-1} { F{A_i^T} F{B_i} }
        Image CoeffA(padded_width, padded_height, 1, 2);
        Image CoeffB(padded_width, padded_height, 1, 2);
        Image Ri[CG_ITERATIONS], Di;
        float alpha = 0.f;
        Image CoeffADi;

        // Compute CoeffB:
        ComplexMultiply::apply(Bx, Px, true);
        ComplexMultiply::apply(By, Py, true);
        ComplexMultiply::apply(dxBx, dxPx, true);
        ComplexMultiply::apply(dyBy, dyPy, true);
        ComplexMultiply::apply(dxyBxy, dxyPxy, true);
        for (int y = 0; y < padded_height; y++)
            for (int x = 0; x < padded_width; x++)
                for (int c = 0; c < 2; c++)
                    CoeffB(x, y)[c] = 25.f * (Bx(x, y)[c] + By(x, y)[c]) + 12.5f * (dxBx(x, y)[c] + dyBy(x, y)[c])
                        + 6.25f * (dxyBxy(x, y)[c]);
        Scale::apply(CoeffB, 2.0f);
        InverseFourierTransform(CoeffB);
    
        // Compute CoeffA
        ComplexMultiply::apply(Px, Px, true);
        ComplexMultiply::apply(Py, Py, true);
        ComplexMultiply::apply(dxPx, dxPx, true);
        ComplexMultiply::apply(dyPy, dyPy, true);
        ComplexMultiply::apply(dxyPxy, dxyPxy, true);
        for (int y = 0; y < padded_height; y++)
            for (int x = 0; x < padded_width; x++)
                for (int c = 0; c < 1; c++) // c=1 should always be zero.
                    CoeffA(x, y)[c] = 25.f * (Px(x, y)[c] + Py(x, y)[c]) + 12.5f * (dxPx(x, y)[c] + dyPy(x, y)[c])
                        + 6.25f * (dxyPxy(x, y)[c]) + beta;
        Scale::apply(CoeffA, 2.0f);

        // Enlarge the kernel
        K = EnlargeKernel(K, padded_width, padded_height);
        toc = currentTime(); printf(" CG Setup  : %.3f sec\n", toc - tic); tic = toc;
    
        // Actual conjugate gradient iterations.
        for (int i = 0; i < CG_ITERATIONS && i < m * m; i++) {
            /*
            sprintf(filename_c, "kernel%d_%d.tmp", iteration, i);
            FileTMP::save(ContractKernel(K, m), std::string(filename_c), "float");
            */

            /*      // Compute f and print it out = 1/2 K^T CoeffA K - CoeffB K.
                    Image K2 = K.copy();
                    FourierTransform(K2);  // K2 = F{K}
                    ComplexMultiply::apply(K2, CoeffA, false); // K2 = F{CoeffA}.*F{K}
                    InverseFourierTransform(K2); // K2 = CoeffA K
                    float score = 0.f;
                    for (int y = 0; y < m; y++) {
                    int y_old = (y - (m >> 1) + padded_height) % padded_height;
                    int x_old = ( - (m >> 1) + padded_width) % padded_width;
                    for (int x = 0; x < m; x++, x_old++) {
                    if (x_old >= padded_width) x_old = 0;
                    score += (K2(x_old, y_old)[0] * 0.5f - CoeffB(x_old, y_old)[0]) * K(x_old, y_old)[0];
                    }
                    }     
                    printf("Iteration %d CG %d: %f\n", iteration, i, score); */

            // 1) Compute residual Ri: CoeffB - CoeffA * K,
            // In subsequent iterations, Ri = CoeffB - CoeffA * (K{i-1} + delta) = R{i-1} - CoeffA * Di * coeff
            Ri[i] = Image(padded_width, padded_height, 1, 1);
            if (i == 0) {
                Image tmp = K.copy();
                FourierTransform(tmp);  // tmp = F{K}
                ComplexMultiply::apply(tmp, CoeffA, false); // tmp =  F{CoeffA} F{K}
                InverseFourierTransform(tmp); // tmp = CoeffA * K
                for (int y = 0; y < m; y++) {
                    int y_old = (y - (m >> 1) + padded_height) % padded_height;
                    int x_old = ( - (m >> 1) + padded_width) % padded_width;
                    for (int x = 0; x < m; x++, x_old++) {
                        if (x_old >= padded_width) x_old = 0;
                        Ri[i](x_old, y_old)[0] = CoeffB(x_old, y_old)[0] - tmp(x_old, y_old)[0];
                    }
                }
            } else {
                for (int y = 0; y < m; y++) {
                    int y_old = (y - (m >> 1) + padded_height) % padded_height;
                    int x_old = ( - (m >> 1) + padded_width) % padded_width;
                    for (int x = 0; x < m; x++, x_old++) {
                        if (x_old >= padded_width) x_old = 0;
                        Ri[i](x_old, y_old)[0] = Ri[i-1](x_old, y_old)[0] - alpha * CoeffADi(x_old, y_old)[0];
                    }
                }	
            }

            // 2) Compute search direction, Di.
            {
                if (i == 0) {
                    Di = Ri[i].copy();
                } else {
                    float modifier_top = 0.f;
                    float modifier_bottom = 0.f;
                    for (int y = 0; y < m; y++) {
                        int y_old = (y - (m >> 1) + padded_height) % padded_height;
                        int x_old = ( - (m >> 1) + padded_width) % padded_width;
                        for (int x = 0; x < m; x++, x_old++) {
                            if (x_old >= padded_width) x_old = 0;
                            modifier_top += Ri[i](x_old, y_old)[0] * Ri[i](x_old, y_old)[0];
                            modifier_bottom += Ri[i-1](x_old, y_old)[0] * Ri[i-1](x_old, y_old)[0];
                        }
                    }
                    float modifier = modifier_top / modifier_bottom;
                    if (modifier_top < 0.000001f) break;
                    for (int y = 0; y < m; y++) {
                        int y_old = (y - (m >> 1) + padded_height) % padded_height;
                        int x_old = ( - (m >> 1) + padded_width) % padded_width;
                        for (int x = 0; x < m; x++, x_old++) {
                            if (x_old >= padded_width) x_old = 0;
                            Di(x_old, y_old)[0] = Ri[i](x_old, y_old)[0] + Di(x_old, y_old)[0] * modifier;
                        }
                    }
                }
            }

            // 3) Compute the coefficient of movement along the search direction Di.
            // It is given by coeff = di^T Ri[0] / di^T CoeffA di
            // Wiki says the numerator should be Ri[i]^T Ri[i]. Those are empirically equal.
            {
                CoeffADi = RealComplex::apply(Di); //  = Di
                FourierTransform(CoeffADi); //  = F{Di}
                ComplexMultiply::apply(CoeffADi, CoeffA, false); //  = F{CoeffA} .* F{Di}
                InverseFourierTransform(CoeffADi); //  = CoeffA Di
                float numerator = 0.f;
                float denominator = 0.f;
                for (int y = 0; y < m; y++) {
                    int y_old = (y - (m >> 1) + padded_height) % padded_height;
                    int x_old = ( - (m >> 1) + padded_width) % padded_width;
                    for (int x = 0; x < m; x++, x_old++) {
                        if (x_old >= padded_width) x_old = 0;
                        numerator += Ri[i](x_old, y_old)[0] * Ri[i](x_old, y_old)[0];
                        denominator += CoeffADi(x_old, y_old)[0] * Di(x_old, y_old)[0];
                    }
                }
                alpha = numerator / denominator;
                for (int y = 0; y < m; y++) {
                    int y_old = (y - (m >> 1) + padded_height) % padded_height;
                    int x_old = ( - (m >> 1) + padded_width) % padded_width;
                    for (int x = 0; x < m; x++, x_old++) {	  
                        if (x_old >= padded_width) x_old = 0;
                        K(x_old, y_old)[0] += Di(x_old, y_old)[0] * alpha;
                    }
                }
            }
        }
        K = ContractKernel(K, m);
        // Threshold K.
        float max_K = 0.f;
        for (int y = 0; y < m; y++)
            for (int x = 0; x < m; x++)
                if (max_K < K(x,y)[0]) max_K = K(x,y)[0];
        for (int y = 0; y < m; y++)
            for (int x = 0; x < m; x++)
                if (max_K * 0.05f > K(x,y)[0]) K(x,y)[0] = 0.f;
        NormalizeSum(K);
        /*    char filename_c[20];
              sprintf(filename_c, "kernel%d_final.tmp", iteration);
              std::string filename(filename_c);
              FileTMP::save(K, filename, "float");*/
        // Recenter the center of mass.
        float avg_x = 0, avg_y = 0;
        for (int y = 0; y < m; y++)
            for (int x = 0; x < m; x++) {
                avg_x += x * K(x,y)[0];
                avg_y += y * K(x,y)[0];
            }
        int offset_x = (int)((avg_x + 0.5f) - (m >> 1));
        int offset_y = (int)((avg_y + 0.5f) - (m >> 1));
        K = Crop::apply(K, offset_x, offset_y, m, m);
        NormalizeSum(K);
        toc = currentTime(); printf(" CG        : %.3f sec\n", toc - tic); tic = toc;

        /**************************************************************/
        /* DECONVOLUTION                                              */
        /**************************************************************/
        if (iteration == kernel_scale.size()) break;
        guess = Deconvolve::applyCho2009(Blurry, K);

        /*
        sprintf(filename_c, "guess%d.tmp", iteration);
        FileTMP::save(guess, std::string(filename_c), "float");
        */

        toc = currentTime(); printf(" Deconvolve: %.3f sec\n", toc - tic); tic = toc;
    }
    return K;
}

#include "footer.h"
#endif
