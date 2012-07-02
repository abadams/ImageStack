#ifndef WLS_H
#define WLS_H
#include "header.h"

class WLS : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, float alpha, float lambda, float tolerance);
};

#include "footer.h"
#endif
