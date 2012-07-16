#ifndef IMAGESTACK_LOCAL_LAPLACIAN_H
#define IMAGESTACK_LOCAL_LAPLACIAN_H
#include "header.h"

class LocalLaplacian : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, float alpha, float beta);
};


#include "footer.h"
#endif
