#ifndef NO_FFTW
#ifndef KERNELESTIMATION_H
#define KERNELESTIMATION_H
#include "header.h"

class KernelEstimation : public Operation {
  public:
    void help();
    void parse(vector<string> args);
    static Image apply(Window im, int kernel_size = 9);

    // Helpers
    static void NormalizeSum(Window im);
    static Image EnlargeKernel(Window im, int w, int h);
    static Image ContractKernel(Window im, int size);
    static Image BilinearResample(Image im, int w, int h);

  private:
    static void ShockFilterIteration(Window im, float dt = 1.f);
    static void BilateralFilterIteration(Window im, float sigma_r);
    static float DotProduct(Window im1, Window im2,
                            int channel = 0, int frame = 0);
};

#include "footer.h"
#endif
#endif
