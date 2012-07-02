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
#include "GaussTransform.h"
#include "Filter.h"
#include "Statistics.h"
#include "header.h"

#define FourierTransform(X) (FFT::apply(X, true, true, false))
#define InverseFourierTransform(X) (IFFT::apply(X, true, true, false))

void KernelEstimation::help() {
    pprintf("-kernelestimation will compute an estimated blur kernel from a blurry"
            " image using the algorithm described in Cho and Lee, 2009. It takes an"
            " optional argument that specifies the kernel size. The default is 25.\n"
            "Usage: ImageStack -load blurred -kernelestimation 25 -deconvolve cho\n");
}

bool KernelEstimation::test() {
    Image dog = Downsample::apply(Load::apply("pics/dog1.jpg"), 2, 2, 1);
    Image kernel(15, 15, 1, 1);
    // Make it a diagonal line
    for (int i = 4; i < 11; i++) {
        kernel(i, i) = 1;
    }
    // With a blob at the center
    kernel(7, 7) = 2;
    FastBlur::apply(kernel, 1, 1, 0);
    normalizeSum(kernel);

    Image blurry = Convolve::apply(dog, kernel);
    Image estimate = KernelEstimation::apply(blurry, kernel.width);

    return nearlyEqual(estimate*20, kernel*20);

}

void KernelEstimation::parse(vector<string> args) {
    assert(args.size() <= 1, "-kernelestimation takes at most one argument.\n");
    int kernelSize = args.size() == 1 ? readInt(args[0]) : 25;
    Image im = apply(stack(0), kernelSize);
    push(im);
}

/*
 * Normalizes the window to have a sum of 1, provided that
 * the window has one channel and one frame, and a nonzero sum.
 */
void KernelEstimation::normalizeSum(Image im) {
    assert(im.channels == 1 && im.frames == 1,
           "The image to be normalized must have one channel and one frame!\n");
    double sum = 0.0;
    for (int y = 0; y < im.height; y++) {
        for (int x = 0; x < im.width; x++) {
            sum += im(x, y);
        }
    }
    im /= sum;
}

Image KernelEstimation::enlargeKernel(Image im, int w, int h) {
    Image ret(w, h, 1, 2); // Add a second complex channel
    for (int y = 0; y < im.height; y++) {
        int newY = (y - (im.height/2) + h) % h;
        int newX =  -(im.width/2) + w;
        for (int x = 0; x < im.width; x++, newX++) {
            if (newX == w) newX = 0;
            ret(newX, newY) = im(x, y);
        }
    }
    return ret;
}

Image KernelEstimation::contractKernel(Image im, int size) {
    Image ret(size, size, 1, 1);
    for (int y = 0; y < size; y++) {
        int yOld = (y - (size / 2) + im.height) % im.height;
        int xOld = (- (size / 2) + im.width) % im.width;
        for (int x = 0; x < size; x++, xOld++) {
            if (xOld >= im.width) xOld = 0;
            ret(x, y) = im(xOld, yOld);
        }
    }
    return ret;
}

/*
 * Applies the shock filter of Osher and Rudin 1990, as described by Cho and Lee 2009.
 */
void KernelEstimation::shockFilterIteration(Image im, float dt) {
    Image input = im.copy();
    for (int t = 0; t < im.frames; t++) {
        for (int y = 1; y < im.height - 1; y++) {
            for (int x = 1; x < im.width - 1; x++) {
                for (int c = 0; c < im.channels; c++) {
                    float dx = (input(x+1, y, t, c) - input(x-1, y, t, c))/2;
                    float dy = (input(x, y+1, t, c) - input(x, y-1, t, c))/2;
                    float laplacian = input(x, y-1, t, c) + input(x, y+1, t, c) +
                                      input(x-1, y, t, c) + input(x+1, y, t, c) - 4 * input(x, y, t, c);
                    im(x, y, t, c) -= sqrtf(dx*dx+dy*dy) * (laplacian > 0.f ? 1.f : -1.f) * dt;
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
void KernelEstimation::bilateralFilterIteration(Image im, float sigmaR) {
    float sigmaS = 2.0f;
    Image gaussian(5, 5, 1, 1);
    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++) {
            float dx = x-2.0f;
            float dy = y-2.0f;
            gaussian(x, y) = expf(-(dx*dx + dy*dy) / (2.0f*sigmaS*sigmaS));
        }
    }
    Image input = im.copy();

    float rangeSigmaMult = -1 / (2.0f * sigmaR * sigmaR);
    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width; x++) {
                    float totalWeight = 0.0f, pixel = 0.0f;
                    for (int dy = -2; dy <= 2; dy++) {
                        if (y + dy < 0) continue;
                        if (y + dy >= im.height) break;
                        for (int dx = -2; dx <= 2; dx++) {
                            if (x + dx < 0) continue;
                            if (x + dx >= im.width) break;
                            float colorDist = input(x+dx, y+dy, t, c) - input(x, y, t, c);
                            float weight = gaussian(dx+2, dy+2);
                            weight *= fastexp(colorDist * colorDist * rangeSigmaMult);
                            totalWeight += weight;
                            pixel += weight * input(x+dx, y+dy, t, c);
                        }
                    }
                    im(x, y, t, c) = pixel / totalWeight;
                }
            }
        }
    }
}

Image KernelEstimation::bilinearResample(Image im, int w, int h) {
    Image ret(w, h, im.frames, im.channels);
    vector<float> sample(im.channels);
    float xFactor = ((float)im.width-1) / (w-1);
    for (int t = 0; t < im.frames; t++) {
        for (int y = 0; y < h - 1; y++) {
            float yOld = ((float)y) / (h-1) * (im.height-1);
            int yOldInt = (int)yOld;
            float yOldFrac = yOld - yOldInt;
            for (int x = 0; x < w - 1; x++) {
                im.sample2DLinear(x * xFactor, yOld, t, sample);
                for (int c = 0; c < im.channels; c++) {
                    ret(x, y, t, c) = sample[c];
                }
            }
            for (int c = 0; c < im.channels; c++) {
                ret(w-1, y, t, c) = im(im.width-1, yOldInt, t, c) * (1 - yOldFrac)
                                    + im(im.width-1, yOldInt + 1, t, c) * (yOldFrac);
            }
        }
        for (int x = 0; x < w - 1; x++) {
            float xOld = x * xFactor;
            int xOldInt = (int)xOld;
            float xOldFrac = xOld - xOldInt;
            for (int c = 0; c < im.channels; c++) {
                ret(x, h-1, t, c) = im(xOldInt, im.height-1, t, c) * (1 - xOldFrac)
                                    + im(xOldInt+1, im.height-1, t, c) * (xOldFrac);
            }
        }
        for (int c = 0; c < im.channels; c++)
            ret(w-1, h-1, t, c) = im(im.width-1, im.height-1, t, c);
    }
    return ret;
}

Image KernelEstimation::apply(Image B, int kernelSize) {

    /******************************* Parameter check */

    assert((kernelSize % 2) == 1, "The desired kernel size must be odd.\n");

    assert((B.channels == 1 || B.channels == 3) && B.frames == 1,
           "The blurred image must have 1 or 3 channels, and 1 frame.\n");

    /******************************* Initialize Constants */
    float dt = 1.f;
    float sigmaR = 0.5f;
    float gradientThreshold = 0.f;
    const int CG_ITERATIONS = 10;
    int primes[] = {2, 3, 5, 7};

    /******************************* Compute Kernel Sizes */
    std::vector<int> kernelScale;
    {
        int tmp = kernelSize;
        kernelScale.push_back(tmp);
        while (tmp != 3) {
            tmp = ((int)((float)tmp / 1.6f) / 2) * 2 + 1;
            kernelScale.push_back(tmp);
        }
    }

    /******************************* Declare Local Variables */
    Image Bgray = (B.channels == 3) ? ColorConvert::apply(B, "rgb", "y") : B.copy();
    Image Blurry;
    Image guess;
    Image K(3, 3, 1, 1);
    K(1,1) = 0.4;
    K(0,1) = K(2,1) = K(1,0) = K(1,2) = 0.15;

    float tic = currentTime(), toc;
    /******************************* Coarse-to-Fine Loop */
    for (unsigned int iteration = 1; iteration <= kernelScale.size(); iteration++) {
        // Setup.
        int m = kernelScale[kernelScale.size()-iteration];
        int newwidth = ((float)m) / kernelSize * B.width;
        int newheight = ((float)m) / kernelSize * B.height;
        int paddedWidth = 0, paddedHeight = 0;
        for (int i = 0; i < 4; i++) {
            int s;
            for (s = 1; s < newwidth + m - 1; s *= primes[i]);
            if (s < paddedWidth || i == 0) paddedWidth = s;
            for (s = 1; s < newheight + m - 1; s *= primes[i]);
            if (s < paddedHeight || i == 0) paddedHeight = s;
        }
        printf("Iteration %d/%d (%dx%d, %dx%d)\n", iteration,
               (int)kernelScale.size(), newwidth, newheight, m, m);

        // Generate the important images at the current scale.
        Blurry = bilinearResample(Bgray, newwidth, newheight);

        //    if (iteration <= 3)
        //      guess = ColorConvert::apply(Load::apply("DeblurTest/input3.tmp"), "rgb", "y");

        guess = bilinearResample(iteration == 1 ? Bgray : guess, newwidth, newheight);

        K = bilinearResample(K, m, m);
        normalizeSum(K);
        toc = currentTime(); printf(" Setup     : %.3f sec\n", toc - tic); tic = toc;

        /**************************************************************/
        /* PREDICTION                                                 */
        /**************************************************************/
        bilateralFilterIteration(guess, sigmaR);
        shockFilterIteration(guess, dt);

        /*char filename_c[20];
        sprintf(filename_c, "prediction%d.tmp", iteration);
        FileTMP::save(guess, std::string(filename_c), "float");*/

        // Compute gradient and generate an angular histogram.
        if (iteration == 1) {
            // Compute the gradient histogram.
            unsigned int gradCap = min((newwidth-1) * (newheight-1), 100000);
            float probabilityCutoff = gradCap / (newwidth-1) * (newheight-1);
            std::vector<float> grad_hist[4];
            for (int i = 0; i < 4; i++) {
                grad_hist[i].clear();
                grad_hist[i].reserve(gradCap);
            }
            for (int y = 0; y < newheight - 1; y++)
                for (int x = 0; x < newwidth - 1; x++) {
                    float dx = guess(x+1, y) - guess(x, y);
                    float dy = guess(x, y+1) - guess(x, y);
                    float angle = atan2(dy, dx);
                    int bin = (int)(floor(angle / (M_PI * 0.25f)) + 4) % 4;
                    if (grad_hist[bin].size() < gradCap &&
                        (double)rand() / RAND_MAX < probabilityCutoff) {
                        grad_hist[bin].push_back(dx*dx+dy*dy);
                    }
                }
            // Find the cutoff.
            float cutoff = 0.f;
            float inclusionCount = m * 4.0f; // m * r in the paper, and r=2.
            for (int i = 0; i < 4; i++) {
                std::sort(grad_hist[i].begin(), grad_hist[i].end());
                int cutoffIndex = 0;
                cutoffIndex = grad_hist[i].size() - 1 - inclusionCount;
                if (cutoffIndex >= 0 && cutoff < grad_hist[i][cutoffIndex])
                    cutoff = grad_hist[i][cutoffIndex];
            }
            gradientThreshold = cutoff * 0.1f;
        } else {
            gradientThreshold *= 0.9f * 0.9f;
        }
        // Generate the gradient channels.
        Image Px(paddedWidth, paddedHeight, 1, 2);
        Image Py(paddedWidth, paddedHeight, 1, 2);
        Image Bx(paddedWidth, paddedHeight, 1, 2);
        Image By(paddedWidth, paddedHeight, 1, 2);
        int xoffset = (paddedWidth - newwidth) / 2;
        int yoffset = (paddedHeight - newheight) / 2;
        for (int y = 0; y < newheight - 1; y++) {
            for (int x = 0; x < newwidth - 1; x++) {
                float dx = guess(x+1, y) - guess(x, y);
                float dy = guess(x, y+1) - guess(x, y);
                if (dx*dx+dy*dy > gradientThreshold) {
                    Px(x + xoffset, y + yoffset) = dx; Py(x + xoffset, y + yoffset) = dy;
                }
                Bx(x + xoffset, y + yoffset) = Blurry(x+1, y) - Blurry(x, y);
                By(x + xoffset, y + yoffset) = Blurry(x, y+1) - Blurry(x, y);
            }
        }
        // Update constants for next iteration.
        dt *= 0.9f;
        sigmaR *= 0.9f;
        toc = currentTime(); printf(" Prediction: %.3f sec\n", toc - tic); tic = toc;
        /*
        sprintf(filenameC, "edges%d.tmp", iteration);
        FileTMP::save(Adjoin::apply(Px, Py, 't'), std::string(filenameC), "float");
        */

        /**************************************************************/
        /* KERNEL ESTIMATION                                          */
        /**************************************************************/
        // Build the gradient images.
        float beta = 1.f;
        Image dxPx(paddedWidth, paddedHeight, 1, 2);
        Image dyPy(paddedWidth, paddedHeight, 1, 2);
        Image dxyPxy(paddedWidth, paddedHeight, 1, 2);
        Image dxBx(paddedWidth, paddedHeight, 1, 2);
        Image dyBy(paddedWidth, paddedHeight, 1, 2);
        Image dxyBxy(paddedWidth, paddedHeight, 1, 2);
        for (int y = yoffset; y < newheight - 2 + yoffset; y++) {
            for (int x = xoffset; x < newwidth - 2 + xoffset; x++) {
                dxPx(x, y) = Px(x+1, y) - Px(x, y);
                dyPy(x, y) = Py(x, y+1) - Py(x, y);
                dxyPxy(x, y) = (Px(x, y+1) - Px(x, y) + Py(x+1, y) - Py(x, y)) * 0.5f;
                dxBx(x, y) = Bx(x+1, y) - Bx(x, y);
                dyBy(x, y) = By(x, y+1) - By(x, y);
                dxyBxy(x, y) = (Bx(x, y+1) - Bx(x, y) + By(x+1, y) - By(x, y)) * 0.5f;
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
        Image CoeffA(paddedWidth, paddedHeight, 1, 2);
        Image CoeffB(paddedWidth, paddedHeight, 1, 2);
        Image Ri[CG_ITERATIONS], Di;
        float alpha = 0.f;
        Image CoeffADi;

        // Compute CoeffB:
        ComplexMultiply::apply(Bx, Px, true);
        ComplexMultiply::apply(By, Py, true);
        ComplexMultiply::apply(dxBx, dxPx, true);
        ComplexMultiply::apply(dyBy, dyPy, true);
        ComplexMultiply::apply(dxyBxy, dxyPxy, true);
        for (int c = 0; c < 2; c++)
            for (int y = 0; y < paddedHeight; y++)
                for (int x = 0; x < paddedWidth; x++)
                    CoeffB(x, y, c) = (50.0f * (Bx(x, y, c) + By(x, y, c)) +
                                       25.0f * (dxBx(x, y, c) + dyBy(x, y, c)) +
                                       12.5f * (dxyBxy(x, y, c)));
        InverseFourierTransform(CoeffB);

        // Compute CoeffA
        ComplexMultiply::apply(Px, Px, true);
        ComplexMultiply::apply(Py, Py, true);
        ComplexMultiply::apply(dxPx, dxPx, true);
        ComplexMultiply::apply(dyPy, dyPy, true);
        ComplexMultiply::apply(dxyPxy, dxyPxy, true);
        for (int y = 0; y < paddedHeight; y++)
            for (int x = 0; x < paddedWidth; x++)
                CoeffA(x, y) = (50.0f * (Px(x, y) + Py(x, y)) +
                                25.0f * (dxPx(x, y) + dyPy(x, y)) +
                                12.5f * (dxyPxy(x, y)) +
                                2.0f * beta);

        // Enlarge the kernel
        K = enlargeKernel(K, paddedWidth, paddedHeight);
        toc = currentTime(); printf(" CG Setup  : %.3f sec\n", toc - tic); tic = toc;

        // Actual conjugate gradient iterations.
        for (int i = 0; i < CG_ITERATIONS && i < m * m; i++) {
            /*
            sprintf(filenameC, "kernel%d_%d.tmp", iteration, i);
            FileTMP::save(contractKernel(K, m), std::string(filenameC), "float");
            */

            /*      // Compute f and print it out = 1/2 K^T CoeffA K - CoeffB K.
                    Image K2 = K.copy();
                    FourierTransform(K2);  // K2 = F{K}
                    ComplexMultiply::apply(K2, CoeffA, false); // K2 = F{CoeffA}.*F{K}
                    InverseFourierTransform(K2); // K2 = CoeffA K
                    float score = 0.f;
                    for (int y = 0; y < m; y++) {
                    int yOld = (y - (m / 2) + paddedHeight) % paddedHeight;
                    int xOld = ( - (m / 2) + paddedWidth) % paddedWidth;
                    for (int x = 0; x < m; x++, xOld++) {
                    if (xOld >= paddedWidth) xOld = 0;
                    score += (K2(xOld, yOld) * 0.5f - CoeffB(xOld, yOld)) * K(xOld, yOld);
                    }
                    }
                    printf("Iteration %d CG %d: %f\n", iteration, i, score); */

            // 1) Compute residual Ri: CoeffB - CoeffA * K,
            // In subsequent iterations, Ri = CoeffB - CoeffA * (K{i-1} + delta) = R{i-1} - CoeffA * Di * coeff
            Ri[i] = Image(paddedWidth, paddedHeight, 1, 1);
            if (i == 0) {
                Image tmp = K.copy();
                FourierTransform(tmp);  // tmp = F{K}
                ComplexMultiply::apply(tmp, CoeffA, false); // tmp =  F{CoeffA} F{K}
                InverseFourierTransform(tmp); // tmp = CoeffA * K
                for (int y = 0; y < m; y++) {
                    int yOld = (y - (m / 2) + paddedHeight) % paddedHeight;
                    int xOld = (- (m / 2) + paddedWidth) % paddedWidth;
                    for (int x = 0; x < m; x++, xOld++) {
                        if (xOld >= paddedWidth) xOld = 0;
                        Ri[i](xOld, yOld) = CoeffB(xOld, yOld) - tmp(xOld, yOld);
                    }
                }
            } else {
                for (int y = 0; y < m; y++) {
                    int yOld = (y - (m / 2) + paddedHeight) % paddedHeight;
                    int xOld = (- (m / 2) + paddedWidth) % paddedWidth;
                    for (int x = 0; x < m; x++, xOld++) {
                        if (xOld >= paddedWidth) xOld = 0;
                        Ri[i](xOld, yOld) = Ri[i-1](xOld, yOld) - alpha * CoeffADi(xOld, yOld);
                    }
                }
            }

            // 2) Compute search direction, Di.
            {
                if (i == 0) {
                    Di = Ri[i].copy();
                } else {
                    float modifierTop = 0.f;
                    float modifierBottom = 0.f;
                    for (int y = 0; y < m; y++) {
                        int yOld = (y - (m / 2) + paddedHeight) % paddedHeight;
                        int xOld = (- (m / 2) + paddedWidth) % paddedWidth;
                        for (int x = 0; x < m; x++, xOld++) {
                            if (xOld >= paddedWidth) xOld = 0;
                            modifierTop += Ri[i](xOld, yOld) * Ri[i](xOld, yOld);
                            modifierBottom += Ri[i-1](xOld, yOld) * Ri[i-1](xOld, yOld);
                        }
                    }
                    float modifier = modifierTop / modifierBottom;
                    if (modifierTop < 0.000001f) break;
                    for (int y = 0; y < m; y++) {
                        int yOld = (y - (m / 2) + paddedHeight) % paddedHeight;
                        int xOld = (- (m / 2) + paddedWidth) % paddedWidth;
                        for (int x = 0; x < m; x++, xOld++) {
                            if (xOld >= paddedWidth) xOld = 0;
                            Di(xOld, yOld) = Ri[i](xOld, yOld) + Di(xOld, yOld) * modifier;
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
                    int yOld = (y - (m / 2) + paddedHeight) % paddedHeight;
                    int xOld = (- (m / 2) + paddedWidth) % paddedWidth;
                    for (int x = 0; x < m; x++, xOld++) {
                        if (xOld >= paddedWidth) xOld = 0;
                        numerator += Ri[i](xOld, yOld) * Ri[i](xOld, yOld);
                        denominator += CoeffADi(xOld, yOld) * Di(xOld, yOld);
                    }
                }
                alpha = numerator / denominator;
                for (int y = 0; y < m; y++) {
                    int yOld = (y - (m / 2) + paddedHeight) % paddedHeight;
                    int xOld = (- (m / 2) + paddedWidth) % paddedWidth;
                    for (int x = 0; x < m; x++, xOld++) {
                        if (xOld >= paddedWidth) xOld = 0;
                        K(xOld, yOld) += Di(xOld, yOld) * alpha;
                    }
                }
            }
        }
        K = contractKernel(K, m);
        // Threshold K.
        float max_K = 0.0f;
        for (int y = 0; y < m; y++) {
            for (int x = 0; x < m; x++) {
                if (max_K < K(x, y)) max_K = K(x,y);
            }
        }
        for (int y = 0; y < m; y++) {
            for (int x = 0; x < m; x++) {
                if (max_K * 0.05f > K(x,y)) K(x,y) = 0.f;
            }
        }
        normalizeSum(K);
        /*    char filenameC[20];
              sprintf(filenameC, "kernel%dFinal.tmp", iteration);
              std::string filename(filenameC);
              FileTMP::save(K, filename, "float");*/
        // Recenter the center of mass.
        float avgX = 0, avgY = 0;
        for (int y = 0; y < m; y++) {
            for (int x = 0; x < m; x++) {
                avgX += x * K(x,y);
                avgY += y * K(x,y);
            }
        }
        float offsetX = avgX - (m / 2);
        float offsetY = avgY - (m / 2);
        K = Translate::apply(K, -offsetX, -offsetY, 0);
        normalizeSum(K);
        toc = currentTime(); printf(" CG        : %.3f sec\n", toc - tic); tic = toc;

        /**************************************************************/
        /* DECONVOLUTION                                              */
        /**************************************************************/
        if (iteration == kernelScale.size()) break;
        guess = Deconvolve::applyCho2009(Blurry, K);

        /*
        sprintf(filenameC, "guess%d.tmp", iteration);
        FileTMP::save(guess, std::string(filenameC), "float");
        */

        toc = currentTime(); printf(" Deconvolve: %.3f sec\n", toc - tic); tic = toc;
    }
    return K;
}

#include "footer.h"
#endif
