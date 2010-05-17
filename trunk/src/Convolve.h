#ifndef IMAGESTACK_CONVOLVE_H
#define IMAGESTACK_CONVOLVE_H
#include "header.h"

class Convolve : public Operation {
  public:
    void help();
    void parse(vector<string> args);

    enum ChannelMode {INNER = 0, OUTER, ELEMENTWISE};
    enum BoundaryCondition {ZERO = 0, HOMOGENEOUS, CLAMP, WRAP};

    static Image apply(Window im, Window filter, BoundaryCondition b = ZERO, ChannelMode m = OUTER);

};


class Deconvolve : public Operation {
 public:
    void help();
    void parse(vector<string> args);
    static Image apply(Window im, Window filter, float maxTime = -1, int maxIterations = -1);
    static Image apply(Window im, Window filter, Window initialGuess, float maxTime = -1, int maxIterations = -1);
 private:
    static double dot(Image &a, Image &b);
    static double norm(Image &a);

    // out = alpha * a + b
    static void scaleAdd(Image out, float alpha, Image &a, Image &b);
};

#include "footer.h"
#endif
