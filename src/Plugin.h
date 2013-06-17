#ifndef IMAGESTACK_PLUGIN_H
#define IMAGESTACK_PLUGIN_H

namespace ImageStack {


class Plugin : public Operation {
public:
    void parse(vector<string> args);
    bool test();
    void help();
};

}

#endif
