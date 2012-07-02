#ifndef IMAGESTACK_COLOR_H
#define IMAGESTACK_COLOR_H
#include "header.h"

class ColorMatrix : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, const vector<float> &matrix);
    static Image apply(Image im, const float *matrix, int outChannels);
};

class ColorConvert : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, string from, string to);
    static Image rgb2hsv(Image im);
    static Image hsv2rgb(Image im);
    static Image rgb2y(Image im);
    static Image y2rgb(Image im);
    static Image rgb2yuv(Image im);
    static Image yuv2rgb(Image im);
    static Image rgb2xyz(Image im);
    static Image xyz2rgb(Image im);
    static Image lab2xyz(Image im);
    static Image xyz2lab(Image im);
    static Image rgb2lab(Image im);
    static Image lab2rgb(Image im);

    static Image uyvy2yuv(Image im);
    static Image yuyv2yuv(Image im);

    static Image uyvy2rgb(Image im);
    static Image yuyv2rgb(Image im);

    static Image argb2xyz(Image im);
    static Image xyz2argb(Image im);

    static Image argb2rgb(Image im);
    static Image rgb2argb(Image im);
};

class Demosaic : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image win, int xoff, int yoff, bool awb);
};

#include "footer.h"
#endif
