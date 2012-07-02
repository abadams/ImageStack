#ifndef IMAGESTACK_HDR_H
#define IMAGESTACK_HDR_H
#include "header.h"

class AssembleHDR : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image frames);
    static Image apply(Image frames, const vector<float> &exposures, string gamma="1.0");
private:
    enum CutoffType { Regular, LongestExposure, ShortestExposure};
    static float weightFunc(float maxVal, CutoffType c = Regular);

};

#include "footer.h"
#endif
