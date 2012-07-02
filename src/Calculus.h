#ifndef IMAGESTACK_CALCULUS_H
#define IMAGESTACK_CALCULUS_H
#include "header.h"

class Gradient : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, string dimensions);
    static void apply(Image im, char dimension);
};

class Integrate : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, string dimensions);
    static void apply(Image im, char dimension);
};

class GradMag : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im);
};

class Poisson : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image dx, Image dy, float termination = 0.01);
};

#include "footer.h"
#endif
