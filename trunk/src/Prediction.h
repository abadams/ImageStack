#ifndef IMAGESTACK_PREDICTION_H
#define IMAGESTACK_PREDICTION_H
namespace ImageStack {

class Inpaint : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, Image mask);
};

class SeamlessClone : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image dst, Image src, Image mask);
    static void apply(Image dst, Image src);
};

}
#endif
