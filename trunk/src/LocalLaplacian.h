#ifndef IMAGESTACK_LOCAL_LAPLACIAN_H
#define IMAGESTACK_LOCAL_LAPLACIAN_H
namespace ImageStack {

class LocalLaplacian : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, float alpha, float beta);
};


}
#endif
