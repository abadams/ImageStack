#ifndef IMAGESTACK_LIGHTFIELD_H
#define IMAGESTACK_LIGHTFIELD_H
#include "header.h"

// a LightField is an image which assumes u and v are rolled up into x
// and y, like an image of the lenslets in a plenoptic camera
class LightField {
public:
    LightField(Image im, int uSize_, int vSize_) : image(im), uSize(uSize_), vSize(vSize_) {
        assert(im.width % uSize == 0, "width is not a multiple of lenslet width\n");
        assert(im.height % vSize == 0, "height is not a multiple of lenslet height\n");
        xSize = im.width / uSize;
        ySize = im.height / vSize;
    }

    float &operator()(int x, int y, int u, int v, int c) {
        return image(x*uSize + u, y*vSize + v, c);
    }

    float &operator()(int x, int y, int u, int v, int t, int c) {
        return image(x*uSize + u, y*vSize + v, t, c);
    }

    // quadrilinear 4D sampling (quadriLanczos3 too expensive, 6^4=1296)
    // x,y,u,v follow the same coordinate conventions as
    // operator()
    void sample4D(float x, float y, float u, float v, int t, float *result) {
        int ix[2], iy[2], iu[2], iv[2]; // integer indices
        float wx[2], wy[2], wu[2], wv[2]; // weighting factors

        if ((x < 0 || y < 0 || x > xSize-1 || y > ySize-1)
            || (u < 0 || v < 0 || u > uSize-1 || v > vSize-1)) {
            // out of bounds, so return zero
            for (int c = 0; c < image.channels; c++) {
                result[c] = 0;
            }
            return;
        }

        ix[0] = (int)(floor(x));
        iy[0] = (int)(floor(y));
        iu[0] = (int)(floor(u));
        iv[0] = (int)(floor(v));
        // clamp against bounds
        ix[0] = clamp(ix[0],0,xSize-1);
        iy[0] = clamp(iy[0],0,ySize-1);
        iu[0] = clamp(iu[0],0,uSize-1);
        iv[0] = clamp(iv[0],0,vSize-1);

        ix[1] = ix[0]+1;
        iy[1] = iy[0]+1;
        iu[1] = iu[0]+1;
        iv[1] = iv[0]+1;
        // clamp against bounds
        ix[1] = min(ix[1],xSize-1);
        iy[1] = min(iy[1],ySize-1);
        iu[1] = min(iu[1],uSize-1);
        iv[1] = min(iv[1],vSize-1);

        // calculate the weights for quadrilinear
        wx[1] = x-ix[0];
        wy[1] = y-iy[0];
        wu[1] = u-iu[0];
        wv[1] = v-iv[0];
        wx[0] = 1-wx[1];
        wy[0] = 1-wy[1];
        wu[0] = 1-wu[1];
        wv[0] = 1-wv[1];

        // do the computation
        for (int c = 0; c < image.channels; c++) {
            result[c] = 0;
        }

        for (int i = 0; i < 2; i++) { // go through iu
            for (int j = 0; j < 2; j++) { // go through ix
                for (int k = 0; k < 2; k++) { // go through iv
                    for (int l = 0; l < 2; l++) { // go through iy
                        for (int c = 0; c < image.channels; c++) {
                            result[c] += ((*this)(ix[j],iy[l],iu[i],iv[k],t,c) *
                                          wx[j]*wy[l]*wu[i]*wv[k]);
                        }
                    }
                }
            }
        }
    }

    void sample4D(float x, float y, float u, float v, float *result) {
        sample4D(x,y,u,v,0,result);
    }

    Image image;
    int uSize, vSize;
    int xSize, ySize;
};

class LFFocalStack : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(LightField im, float minAlpha, float maxAlpha, float deltaAlpha);
};

class LFPoint : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(LightField lf, float x, float y, float z);
};

#include "footer.h"
#endif
