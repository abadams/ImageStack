#ifndef IMAGESTACK_DISPLAY_H
#define IMAGESTACK_DISPLAY_H
namespace ImageStack {

class Display : public Operation {
public:
    ~Display();
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, bool fullscreen = false);
};


}
#endif
