#ifndef IMAGESTACK_CONVOLVE_H
#define IMAGESTACK_CONVOLVE_H

#include "Arithmetic.h"

#include "header.h"

class Convolve : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);

    enum BoundaryCondition {Zero = 0, Homogeneous, Clamp, Wrap};

    static Image apply(Image im, Image filter, BoundaryCondition b = Zero,
                       Multiply::Mode m = Multiply::Outer);
private:
    static void convolveSingle(Image im, Image filter, Image out, BoundaryCondition b);
};

#include "footer.h"
#endif
