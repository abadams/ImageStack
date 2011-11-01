#ifndef IMAGESTACK_OPTICALFLOW_H
#define IMAGESTACK_OPTICALFLOW_H
#include "header.h"

class DenseCorrespondence : public Operation {
 public:
  void help();
  void parse(vector<string> args);
  static Image apply(Window source, Window target,
		     const int from=0, const int to=0);
 private:
};

class OpticalFlow : public Operation {
public:
    void help();
    void parse(vector<string> args);
    static Image apply(Window source, Window target, Window initial);
    static Image apply(Window source, Window target);

private:
    static bool IsDisplay;
    static void Coarse2FineFlow(Image &vx, Image &vy, Image &warpI2, const Image Im1, const Image Im2, double alpha, double ratio, int minWidth, int nOuterFPIterations, int nInnerFPIterations, int nCGIterations);
    static void RefineFlow(Image initial, Image &vx, Image &vy, Image &warpI2, const Image Im1, const Image Im2, double alpha, double ratio, int minWidth, int nOuterFPIterations, int nInnerFPIterations, int nCGIterations);
    static Image im2feature(Image im);
    static void SmoothFlowPDE(const Image Im1, const Image Im2, Image &warpIm2, Image &u, Image &v, double alpha, int nOuterFPIterations, int nInnerFPIterations, int nCGIterations);
    static void getDxs(Image &imdx, Image &imdy, Image &imdt, Image im1, Image im2);
    static void genInImageMask(Image &mask, Image vx, Image vy);
    static void reset(Image &im);
    static void warpFL(Image &warpIm2, Image Im1, Image Im2, Image vx, Image vy);
    static Image collapse(Image in);
    static void Laplacian(Image &output, Image input, Image weight);
    static double norm2(Image im);
    static double innerproduct(Image im1, Image im2);
    static Image add(Image a, Image b);
    static Image subtract(Image a, Image b);
    static Image gradient(Image im, char dimension);
    static Image multiply(Image a, Image b, Image c);
    static Image multiply(Image a, Image b);
    static Image addAfterScale(Image a, Image b, double c);
    static Image phi(Image a);
    static Image evalConfidence(Image source, Image target, Image flow, float alpha = 0.02, float gamma = 1.0);
};

class OpticalFlowWarp : public Operation {
public:
    void help();
    void parse(vector<string> args);
    static Image apply(Window input, Window from, Window to);
};

#include "footer.h"
#endif
