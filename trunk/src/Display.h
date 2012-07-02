#ifndef IMAGESTACK_DISPLAY_H
#define IMAGESTACK_DISPLAY_H
#include "header.h"

class Display : public Operation {
public:
    ~Display();
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, bool fullscreen = false);
};


#include "footer.h"
#endif
