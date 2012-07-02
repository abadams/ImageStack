#ifndef IMAGESTACK_GAUSS_TRANSFORM_H
#define IMAGESTACK_GAUSS_TRANSFORM_H
#include "header.h"


class GaussTransform : public Operation {
public:
    enum Method {AUTO = 0, EXACT, GRID, PERMUTOHEDRAL, GKDTREE};

    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image slicePositions, Image splatPositions, Image values,
                       vector<float> sigmas, Method m = AUTO);
};


class JointBilateral : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image image, Image reference,
                      float filterWidth, float filterHeight, float filterFrames, float colorSigma,
                      GaussTransform::Method m = GaussTransform::AUTO);
};

class Bilateral : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image image, float filterWidth, float filterHeight,
                      float filterFrames, float colorSigma,
                      GaussTransform::Method m = GaussTransform::AUTO);
};


class BilateralSharpen : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, float spatialSigma, float colorSigma, float sharpness);
};

class ChromaBlur : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, float spatialSigma, float colorSigma);
};



class NLMeans : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image image, float patchSize, int dimensions,
                      float spatialSigma, float patchSigma,
                      GaussTransform::Method m = GaussTransform::AUTO);
};

class FastNLMeans : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image image, float patchSize, float spatialSigma, float patchSigma);
};

class NLMeans3D : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image image, float patchSize, int dimensions,
                      float spatialSigma, float patchSigma,
                      GaussTransform::Method m = GaussTransform::AUTO);
};

#include "footer.h"
#endif
