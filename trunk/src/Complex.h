#ifndef IMAGESTACK_COMPLEX_H
#define IMAGESTACK_COMPLEX_H
#include "header.h"

class ComplexMultiply : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image a, Image b, bool conj = false);
};

class ComplexDivide : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image a, Image b, bool conj = false);
};

class ComplexReal : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im);
};

class RealComplex : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im);
};

class ComplexImag : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im);
};

class ComplexMagnitude : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im);
};

class ComplexPhase : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im);
};

class ComplexConjugate : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im);
};

#include "footer.h"
#endif
