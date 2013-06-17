#ifndef WLS_H
#define WLS_H
namespace ImageStack {

class WLS : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, float alpha, float lambda, float tolerance);
};

}
#endif
