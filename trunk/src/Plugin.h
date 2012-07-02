#ifndef IMAGESTACK_PLUGIN_H
#define IMAGESTACK_PLUGIN_H

#include "header.h"


class Plugin : public Operation {
public:
    void parse(vector<string> args);
    bool test();
    void help();
};

#include "footer.h"

#endif
