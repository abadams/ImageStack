#ifndef IMAGESTACK_PAINT_H
#define IMAGESTACK_PAINT_H
#include "header.h"

class Eval : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, string expression);
};

class EvalChannels : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, vector<string> expressions);
};

class Plot : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, int width, int height, float lineThickness);
};

class Composite : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image dst, Image src);
    static void apply(Image dst, Image src, Image mask);
};

#include "footer.h"
#endif
