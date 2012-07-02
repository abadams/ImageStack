#ifndef IMAGESTACK_GEOMETRY_H
#define IMAGESTACK_GEOMETRY_H
#include "header.h"

class Upsample : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, int boxWidth, int boxHeight, int boxFrames = 1);
};

class Downsample : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, int boxWidth, int boxHeight, int boxFrames = 1);
};

class Subsample : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, int boxWidth, int boxHeight, int boxFrames,
                       int offsetX, int offsetY, int offsetT);
    static Image apply(Image im, int boxWidth, int boxHeight,
                       int offsetX, int offsetY);
    static Image apply(Image im, int boxFrames, int offsetT);
};

class Interleave : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, int rx, int ry, int rt = 1);
};

class Deinterleave : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, int ix, int iy, int it = 1);
};

class Resample : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, int width, int height);
    static Image apply(Image im, int width, int height, int frames);
private:
    static void computeWeights(int oldSize, int newSize, vector<vector<pair<int, float> > > &matrix);
    static Image resampleT(Image im, int frames);
    static Image resampleX(Image im, int width);
    static Image resampleY(Image im, int height);
};

class Rotate : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, float degrees);
};

class AffineWarp : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, vector<float> warp);
    static Image apply(Image im, float *warp);
};

class Crop : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, int minX, int minY, int width, int height);
    static Image apply(Image im, int minX, int minY, int minT, int width, int height, int frames);
    static Image apply(Image im);
};

class Flip : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, char dimension);
};

class Adjoin : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image a, Image b, char dimension);
};

class Transpose : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, char arg1, char arg2);
};

class Translate : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, float xoff, float yoff, float toff = 0);
private:
    static Image applyX(Image im, float xoff);
    static Image applyY(Image im, float yoff);
    static Image applyT(Image im, float toff);
};

class Paste : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image into, Image from,
                      int xdst, int ydst,
                      int xsrc, int ysrc,
                      int width, int height);

    static void apply(Image into, Image from,
                      int xdst, int ydst, int tdst = 0);

    static void apply(Image into, Image from,
                      int xdst, int ydst, int tdst,
                      int xsrc, int ysrc, int tsrc,
                      int width, int height, int frames);
};

class Tile : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, int xTiles, int yTiles, int tTiles = 1);
};

class TileFrames : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, int xTiles, int yTiles);
};

class FrameTiles : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, int xTiles, int yTiles);
};

class Warp : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image coords, Image source);
};

class Reshape : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, int newWidth, int newHeight, int newFrames, int newChannels);
};

#include "footer.h"
#endif
