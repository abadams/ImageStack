#ifndef IMAGESTACK_MATH_H
#define IMAGESTACK_MATH_H
#include "header.h"

class Add : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
};

class Multiply : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    enum Mode {Elementwise = 0, Inner, Outer};
    static Image apply(Image a, Image b, Mode m);
};

class Subtract : public Operation {
public:
    void help();
    void parse(vector<string> args);
    bool test();
};

class Divide : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
};

class Maximum : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image a, Image b);
};

class Minimum : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image a, Image b);
};

class Log : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image a);
};

class Exp : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image a, float base = E);
};

class Abs : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image a);
};

class Offset : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
};

class Scale : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
};

class Gamma : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image a, float);
};

class Mod : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image a, float);
};

class Clamp : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image a, float lower = 0, float upper = 1);
};

class DeNaN : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image a, float replacement = 0);
};

class Threshold : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image a, float val);
};

class Normalize : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image a);
};

class Quantize : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image a, float increment);
};

#include "footer.h"
#endif
