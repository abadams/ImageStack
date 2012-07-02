#ifndef IMAGESTACK_LOCAL_LAPLACIAN_H
#define IMAGESTACK_LOCAL_LAPLACIAN_H
#include "header.h"

class LocalLaplacian : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, float alpha, float beta);
private:
    static Image pyramidDown(Image im);
    static Image pyramidUp(Image im, int w, int h, int f);
};


#include "footer.h"
#endif
