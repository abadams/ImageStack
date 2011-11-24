#ifndef NO_FFTW
#ifndef DECONVOLUTION_H
#define DECONVOLUTION_H
#include "header.h"

class Deconvolve : public Operation {
  public:
    void help();
    void parse(vector<string> args);
    static Image applyCho2009(Window im, Window kernel);
    static Image applyShan2008(Window im, Window kernel);
    static Image applyLevin2007(Window im, Window kernel, float weight);
  private:
    static Image applyPadding(Window im);
};

#include "footer.h"
#endif
#endif
