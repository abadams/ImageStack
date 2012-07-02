#ifndef NO_FFTW
#ifndef DECONVOLUTION_H
#define DECONVOLUTION_H
#include "header.h"

class Deconvolve : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image applyCho2009(Image im, Image kernel);
    static Image applyShan2008(Image im, Image kernel);
    static Image applyLevin2007(Image im, Image kernel, float weight);
private:
    static Image applyPadding(Image im);
};

#include "footer.h"
#endif
#endif
