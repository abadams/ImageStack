#ifndef IMAGESTACK_WAVELET_H
#define IMAGESTACK_WAVELET_H
#include "header.h"

class Haar : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, int times = -1);
};

class InverseHaar : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, int times = -1);
};

class Daubechies : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im);
};

class InverseDaubechies : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im);
};

#include "footer.h"
#endif
