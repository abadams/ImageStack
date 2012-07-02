#ifndef NO_FFTW
#ifndef KERNELESTIMATION_H
#define KERNELESTIMATION_H
#include "header.h"

class KernelEstimation : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, int kernelSize = 25);

    // Helpers
    static void normalizeSum(Image im);
    static Image enlargeKernel(Image im, int w, int h);
    static Image contractKernel(Image im, int size);
    static Image bilinearResample(Image im, int w, int h);

private:
    static void shockFilterIteration(Image im, float dt = 1.0f);
    static void bilateralFilterIteration(Image im, float sigmaR);
};

#include "footer.h"
#endif
#endif
