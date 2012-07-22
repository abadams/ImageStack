#ifndef IMAGESTACK_IMAGE_H
#define IMAGESTACK_IMAGE_H

#include "Expr.h"

#include "tables.h"
#include "header.h"

// The image data type.

// It's a reference-counted pointer type.

// Note that "const Image" means that the reference doesn't change, not
// that the pixel data doesn't. It's equivalent to "float * const
// foo", not "const float * foo". This means that methods tagged const
// are those which do not change the metadata, not those which do not
// change the pixel data.


template<typename SX, typename SY, typename ST, typename SC, bool AffineCase, bool ShiftedCase>
class ImageRef;

class Image {
public:

    int width, height, frames, channels;
    int ystride, tstride, cstride;

    Image() :
        width(0), height(0), frames(0), channels(0),
        ystride(0), tstride(0), cstride(0), data(), base(NULL) {
    }

    Image(int w, int h, int f, int c) :
        width(w), height(h), frames(f), channels(c),
        ystride(w), tstride(w * h), cstride(w * h * f),
        data(new Payload(w * h * f * c + 16)), base(compute_base(data)) {
        // + 16 so that we can walk forwards up to one avx vector
        // width (8 floats) for alignment, and so that we have at
        // least one vector worth of allocation beyond the end so that
        // we can pull vectors safely from the image even if they go
        // off the end
    }

    inline float &operator()(int x) const {
        return (*this)(x, 0, 0, 0);
    }

    inline float &operator()(int x, int y) const {
        return (*this)(x, y, 0, 0);
    }

    inline float &operator()(int x, int y, int c) const {
        return (*this)(x, y, 0, c);
    }

    inline float &operator()(int x, int y, int t, int c) const {
#ifdef BOUNDS_CHECKING
        assert(x >= 0 && x < width &&
               y >= 0 && y < height &&
               t >= 0 && t < frames &&
               c >= 0 && c < channels,
               "Access out of bounds: %d %d %d %d\n",
               x, y, t, c);
#endif
        return (((base + c*cstride) + t*tstride) + y*ystride)[x];
    }

    float *baseAddress() const {
        return base;
    }

    Image copy() const {
        Image m(width, height, frames, channels);
        m.set(*this);
        return m;
    }

    const Image region(int x, int y, int t, int c,
                       int xs, int ys, int ts, int cs) const {
        return Image(*this, x, y, t, c, xs, ys, ts, cs);
    }

    const Image column(int x) const {
        return region(x, 0, 0, 0, 1, height, frames, channels);
    }

    const Image row(int y) const {
        return region(0, y, 0, 0, width, 1, frames, channels);
    }

    const Image frame(int t) const {
        return region(0, 0, t, 0, width, height, 1, channels);
    }

    const Image channel(int c) const {
        return region(0, 0, 0, c, width, height, frames, 1);
    }

    const Image selectColumns(int x, int s) {
        return region(x, 0, 0, 0, s, height, frames, channels);
    }

    const Image selectRows(int x, int s) {
        return region(0, x, 0, 0, width, s, frames, channels);
    }

    const Image selectFrames(int x, int s) {
        return region(0, 0, x, 0, width, height, s, channels);
    }

    const Image selectChannels(int x, int s) {
        return region(0, 0, 0, x, width, height, frames, s);
    }

    bool dense() const {
        return (cstride == width *height *frames && tstride == width *height && ystride == width);
    }


    bool defined() const {
        return base != NULL;
    }

    bool operator==(const Image &other) const {
        return (base == other.base &&
                ystride == other.ystride &&
                tstride == other.tstride &&
                cstride == other.cstride &&
                width == other.width &&
                height == other.height &&
                frames == other.frames &&
                channels == other.channels);
    }

    bool operator!=(const Image &other) const {
        return !(*this == other);
    }

    template<typename A, typename B, typename Enable = FloatExprType(A) >
    struct ExprCheck {
        typedef B result;
    };

    template<typename T>
    typename ExprCheck<T, void>::result operator+=(const T &other) const {
        set((*this) + other);
    }

    template<typename T>
    typename ExprCheck<T, void>::result operator*=(const T &other) const {
        set((*this) * other);
    }

    template<typename T>
    typename ExprCheck<T, void>::result operator-=(const T &other) const {
        set((*this) - other);
    }

    template<typename T>
    typename ExprCheck<T, void>::result operator/=(const T &other) const {
        set((*this) / other);
    }

    typedef enum {ZERO = 0, NEUMANN} BoundaryCondition;

    void sample2D(float fx, float fy, int t, vector<float> &result, BoundaryCondition boundary = ZERO) const {
        sample2D(fx, fy, t, &result[0], boundary);
    }

    void sample2D(float fx, float fy, int t, float *result, BoundaryCondition boundary = ZERO) const {
        int ix = (int)fx;
        int iy = (int)fy;
        const int LEFT = -2;
        const int RIGHT = 3;
        const int WIDTH = 6;
        int minX = ix + LEFT;
        int maxX = ix + RIGHT;
        int minY = iy + LEFT;
        int maxY = iy + RIGHT;

        float weightX[WIDTH];
        float weightY[WIDTH];
        float totalXWeight = 0, totalYWeight = 0;
        for (int x = 0; x < WIDTH; x++) {
            float diff = (fx - (x + ix + LEFT)); // ranges between +/- RIGHT
            float val = lanczos_3(diff);
            weightX[x] = val;
            totalXWeight += val;
        }

        for (int y = 0; y < WIDTH; y++) {
            float diff = (fy - (y + iy + LEFT)); // ranges between +/- RIGHT
            float val = lanczos_3(diff);
            weightY[y] = val;
            totalYWeight += val;
        }

        totalXWeight = 1.0f/totalXWeight;
        totalYWeight = 1.0f/totalYWeight;

        for (int i = 0; i < WIDTH; i++) {
            weightX[i] *= totalXWeight;
            weightY[i] *= totalYWeight;
        }

        for (int c = 0; c < channels; c++) {
            result[c] = 0;
        }

        if (boundary == NEUMANN) {

            float *yWeightPtr = weightY;
            for (int y = minY; y <= maxY; y++) {
                float *xWeightPtr = weightX;
                int sampleY = clamp(0, y, height-1);
                for (int x = minX; x <= maxX; x++) {
                    int sampleX = clamp(0, x, width-1);
                    float yxWeight = (*yWeightPtr) * (*xWeightPtr);
                    for (int c = 0; c < channels; c++) {
                        result[c] += (*this)(sampleX, sampleY, t, c) * yxWeight;
                    }
                    xWeightPtr++;
                }
                yWeightPtr++;
            }
        } else {
            float *weightYBase = weightY;
            float *weightXBase = weightX;
            if (minY < 0) {
                weightYBase -= minY;
                minY = 0;
            }
            if (minX < 0) {
                weightXBase -= minX;
                minX = 0;
            }
            if (maxX > width-1) { maxX = width-1; }
            if (maxY > height-1) { maxY = height-1; }
            float *yWeightPtr = weightYBase;
            for (int y = minY; y <= maxY; y++) {
                float *xWeightPtr = weightXBase;
                for (int x = minX; x <= maxX; x++) {
                    float yxWeight = (*yWeightPtr) * (*xWeightPtr);
                    for (int c = 0; c < channels; c++) {
                        result[c] += (*this)(x, y, t, c) * yxWeight;
                    }
                    xWeightPtr++;
                }
                yWeightPtr++;
            }

        }
    }

    void sample2D(float fx, float fy, vector<float> &result) const {
        sample2D(fx, fy, 0, result);
    }

    void sample2D(float fx, float fy, float *result) const {
        sample2D(fx, fy, 0, result);
    }

    void sample2DLinear(float fx, float fy, vector<float> &result) const {
        sample2DLinear(fx, fy, 0, result);
    }

    void sample2DLinear(float fx, float fy, float *result) const {
        sample2DLinear(fx, fy, 0, result);
    }

    void sample2DLinear(float fx, float fy, int t, vector<float> &result) const {
        sample2DLinear(fx, fy, t, &result[0]);
    }

    void sample2DLinear(float fx, float fy, int t, float *result) const {
        int ix = (int)fx;
        int iy = (int)fy;
        fx -= ix;
        fy -= iy;

        for (int c = 0; c < channels; c++) {
            float s1 = (1-fx) * (*this)(ix, iy, t, c) + fx * (*this)(ix+1, iy, t, c);
            float s2 = (1-fx) * (*this)(ix, iy+1, t, c) + fx * (*this)(ix+1, iy+1, t, c);
            result[c] = (1-fy) * s1 + fy * s2;
        }

    }

    void sample3DLinear(float fx, float fy, float ft, vector<float> &result) const {
        sample3DLinear(fx, fy, ft, &result[0]);
    }

    void sample3DLinear(float fx, float fy, float ft, float *result) const {
        int ix = (int)fx;
        int iy = (int)fy;
        int it = (int)ft;
        fx -= ix;
        fy -= iy;
        ft -= it;

        for (int c = 0; c < channels; c++) {
            float s11 = (1-fx) * (*this)(ix, iy, it, c) + fx * (*this)(ix+1, iy, it, c);
            float s12 = (1-fx) * (*this)(ix, iy+1, it, c) + fx * (*this)(ix+1, iy+1, it, c);
            float s1 = (1-fy) * s11 + fy * s12;

            float s21 = (1-fx) * (*this)(ix, iy, it+1, c) + fx * (*this)(ix+1, iy, it+1, c);
            float s22 = (1-fx) * (*this)(ix, iy+1, it+1, c) + fx * (*this)(ix+1, iy+1, it+1, c);
            float s2 = (1-fy) * s21 + fy * s22;

            result[c] = (1-ft) * s1 + ft * s2;
        }

    }

    void sample3D(float fx, float fy, float ft,
                  vector<float> &result, BoundaryCondition boundary = ZERO) const {
        sample3D(fx, fy, ft, &result[0], boundary);
    }

    void sample3D(float fx, float fy, float ft,
                  float *result, BoundaryCondition boundary = ZERO) const {
        int ix = (int)fx;
        int iy = (int)fy;
        int it = (int)ft;
        const int LEFT = -2;
        const int RIGHT = 3;
        const int WIDTH = 6;
        int minX = ix + LEFT;
        int maxX = ix + RIGHT;
        int minY = iy + LEFT;
        int maxY = iy + RIGHT;
        int minT = it + LEFT;
        int maxT = it + RIGHT;
        float weightX[WIDTH];
        float weightY[WIDTH];
        float weightT[WIDTH];

        float totalXWeight = 0, totalYWeight = 0, totalTWeight = 0;

        for (int x = 0; x < WIDTH; x++) {
            float diff = (fx - (x + ix + LEFT)); // ranges between +/- RIGHT
            float val = lanczos_3(diff);
            weightX[x] = val;
            totalXWeight += val;
        }

        for (int y = 0; y < WIDTH; y++) {
            float diff = (fy - (y + iy + LEFT)); // ranges between +/- RIGHT
            float val = lanczos_3(diff);
            weightY[y] = val;
            totalYWeight += val;
        }

        for (int t = 0; t < WIDTH; t++) {
            float diff = (ft - (t + it + LEFT)); // ranges between +/- RIGHT
            float val = lanczos_3(diff);
            weightT[t] = val;
            totalTWeight += val;
        }

        totalXWeight = 1.0f/totalXWeight;
        totalYWeight = 1.0f/totalYWeight;
        totalTWeight = 1.0f/totalTWeight;

        for (int i = 0; i < WIDTH; i++) {
            weightX[i] *= totalXWeight;
            weightY[i] *= totalYWeight;
            weightT[i] *= totalTWeight;
        }

        for (int c = 0; c < channels; c++) {
            result[c] = 0;
        }

        if (boundary == NEUMANN) {

            float *tWeightPtr = weightT;
            for (int t = minT; t <= maxT; t++) {
                int sampleT = clamp(t, 0, frames-1);
                float *yWeightPtr = weightY;
                for (int y = minY; y <= maxY; y++) {
                    int sampleY = clamp(y, 0, height-1);
                    float tyWeight = (*yWeightPtr) * (*tWeightPtr);
                    float *xWeightPtr = weightX;
                    for (int x = minX; x <= maxX; x++) {
                        int sampleX = clamp(x, 0, width-1);
                        float tyxWeight = tyWeight * (*xWeightPtr);
                        for (int c = 0; c < channels; c++) {
                            result[c] += (*this)(sampleX, sampleY, sampleT, c) * tyxWeight;
                        }
                        xWeightPtr++;
                    }
                    yWeightPtr++;
                }
                tWeightPtr++;
            }

        } else {

            float *weightTBase = weightT;
            float *weightYBase = weightY;
            float *weightXBase = weightX;

            if (minY < 0) {
                weightYBase -= minY;
                minY = 0;
            }
            if (minX < 0) {
                weightXBase -= minX;
                minX = 0;
            }
            if (minT < 0) {
                weightTBase -= minT;
                minT = 0;
            }
            if (maxX > width-1) { maxX = width-1; }
            if (maxY > height-1) { maxY = height-1; }
            if (maxT > frames-1) { maxT = frames-1; }

            float *tWeightPtr = weightTBase;
            for (int t = minT; t <= maxT; t++) {
                float *yWeightPtr = weightYBase;
                for (int y = minY; y <= maxY; y++) {
                    float *xWeightPtr = weightXBase;
                    for (int x = minX; x <= maxX; x++) {
                        float yxWeight = (*yWeightPtr) * (*xWeightPtr);
                        for (int c = 0; c < channels; c++) {
                            result[c] += (*this)(x, y, t, c) * yxWeight;
                        }
                        xWeightPtr++;
                    }
                    yWeightPtr++;
                }
                tWeightPtr++;
            }
        }

    }

    // Evaluate a expression object defined in Expr.h
    // The second argument prevents calls to set from things that are
    // neither floats nor float expression types
    template<typename T>
    void set(const T expr_, const FloatExprType(T) *check = NULL) const {
        FloatExprType(T) expr(expr_);
        {
            assert(defined(), "Can't set undefined image\n");
            int w = expr.getSize(0), h = expr.getSize(1),
                f = expr.getSize(2), c = expr.getSize(3);
            assert((w == 0 || width == w) &&
                   (h == 0 || height == h) &&
                   (f == 0 || frames == f) &&
                   (c == 0 || channels == c),
                   "Can only assign from source of matching size\n");
        }
        
        // Compute the domain over which we can vectorize
        bool boundedVX = expr.boundedVecX();
        int minVX = expr.minVecX();
        int maxVX = expr.maxVecX();


        
        Expr::Region r = {0, 0, 0, 0, width, height, frames, channels};

        // Figure out what regions of what functions are required here
        //float t1 = currentTime();
        expr.prepare(r, 0);
        //float t2 = currentTime();
        expr.prepare(r, 1);
        //float t3 = currentTime();
        expr.prepare(r, 2);
        //float t4 = currentTime();
        

        for (int c = 0; c < channels; c++) {
            for (int t = 0; t < frames; t++) {
                #ifdef _OPENMP
                #pragma omp parallel for
                #endif
                for (int y = 0; y < height; y++) {
                    //printf("Evaluating at scanline %d\n", y);
                    FloatExprType(T)::Iter iter = expr.scanline(0, y, t, c, width);
                    float *const dst = base + c*cstride + t*tstride + y*ystride;
                    ImageStack::Expr::setScanline(iter, dst, 0, width, boundedVX, minVX, maxVX);
                }
            }
        }
        //float t5 = currentTime();

        // Clean up any resources
        expr.prepare(r, 3);
        //float t6 = currentTime();
        //printf("%f %f %f %f %f\n", t2-t1, t3-t2, t4-t3, t5-t4, t6-t5);
    }

    void set(const Expr::Func &func) const {
        realizeFuncIntoImage(*this, func);
    }
    Image(const Expr::Func &func) {
        (*this) = realizeFuncIntoNewImage(func);
    }


    // A version of set that takes a set of up to 4 expressions and
    // sets the image's channels accordingly. This is more efficient
    // than calling set repeatedly when the expressions share common
    // subexpressions. At each pixel, each expression is evaluated,
    // and then each value is stored, preventing nasty chicken-and-egg
    // problems when, for example, permuting channels.
    template<typename A, typename B, typename C, typename D>
    void setChannels(const A a, const B b, const C c, const D d,
                     FloatExprType(A) *pa = NULL,
                     FloatExprType(B) *pb = NULL,
                     FloatExprType(C) *pc = NULL,
                     FloatExprType(D) *pd = NULL) const {
        setChannelsGeneric<4, FloatExprType(A), FloatExprType(B), FloatExprType(C), FloatExprType(D)>(
            FloatExprType(A)(a),
            FloatExprType(B)(b),
            FloatExprType(C)(c),
            FloatExprType(D)(d));
    }

    template<typename A, typename B, typename C>
    void setChannels(const A &a, const B &b, const C &c,
                     FloatExprType(A) *pa = NULL,
                     FloatExprType(B) *pb = NULL,
                     FloatExprType(C) *pc = NULL) const {
        setChannelsGeneric<3, FloatExprType(A), FloatExprType(B), FloatExprType(C), FloatExprType(float)>(
            FloatExprType(A)(a),
            FloatExprType(B)(b),
            FloatExprType(C)(c),
            FloatExprType(float)(0));
    }

    template<typename A, typename B>
    void setChannels(const A &a, const B &b,
                     FloatExprType(A) *pa = NULL,
                     FloatExprType(B) *pb = NULL) const {
        setChannelsGeneric<2, FloatExprType(A), FloatExprType(B), FloatExprType(float), FloatExprType(float)>(
            FloatExprType(A)(a),
            FloatExprType(B)(b),
            FloatExprType(float)(0),
            FloatExprType(float)(0));
    }

#define ShiftedCase(SX, SY, ST, SC)                                     \
    (ShiftedInX(SX) && IndependentOfX(SY) && IndependentOfX(ST) && IndependentOfX(SC)) 
        
#define AffineCase(SX, SY, ST, SC)                                      \
    (AffineInX(SX) && IndependentOfX(SY) && IndependentOfX(ST) && IndependentOfX(SC)) 


    template<typename SX, typename SY, typename ST, typename SC>
    ImageRef<IntExprType(SX), IntExprType(SY), IntExprType(ST), IntExprType(SC), 
             AffineCase(SX, SY, ST, SC),
             ShiftedCase(SX, SY, ST, SC)>
    operator()(const SX &x, const SY &y, const ST &t, const SC &c) const;

    template<typename SX, typename SY, typename SC>
    ImageRef<IntExprType(SX), IntExprType(SY), Expr::ConstInt, IntExprType(SC), 
             AffineCase(SX, SY, Expr::ConstInt, SC), 
             ShiftedCase(SX, SY, Expr::ConstInt, SC)>
    operator()(const SX &x, const SY &y, const SC &c) const;
    
    template<typename SX, typename SY>
    ImageRef<IntExprType(SX), IntExprType(SY), Expr::ConstInt, Expr::ConstInt, 
             AffineCase(SX, SY, Expr::ConstInt, Expr::ConstInt),
             ShiftedCase(SX, SY, Expr::ConstInt, Expr::ConstInt)>
    operator()(const SX &x, const SY &y) const;

    template<typename SX>
    ImageRef<IntExprType(SX), Expr::ConstInt, Expr::ConstInt, Expr::ConstInt, 
             AffineCase(SX, Expr::ConstInt, Expr::ConstInt, Expr::ConstInt),
             ShiftedCase(SX, Expr::ConstInt, Expr::ConstInt, Expr::ConstInt)>
    operator()(const SX &x) const;

    // An image itself is one such expression thing. Here's the
    // interface it needs to implement for that to happen:
    typedef Image FloatExpr;
    const static bool dependsOnX = true;
    int getSize(int i) const {
        switch (i) {
        case 0: return width;
        case 1: return height;
        case 2: return frames;
        case 3: return channels;
        }
        return 0;
    }
    #ifdef BOUNDS_CHECKING
    struct Iter {
        int width;
        const float *const addr;
        Iter() : addr(NULL) {}
        Iter(const float *a, int w) : width(w), addr(a) {}
        float operator[](int x) const {
            assert(x >= 0 && x < width, 
                   "Access out of bounds in image iterator:\n"
                   "%d is not within 0 - %d\n", x, width);
            return addr[x];
        }
        ImageStack::Expr::Vec::type vec(int x) const {
            assert(x >= 0 && x <= width - ImageStack::Expr::Vec::width,
                   "Vector access out of bounds in image iterator:\n"
                   "%d is not sufficiently within 0 - %d\n", x, width);
            return ImageStack::Expr::Vec::load(addr+x);
        }
    };
    Iter scanline(int x, int y, int t, int c, int w) const {
        assert(x >= 0 && x+w <= width,
               "Scanline will access image out of bounds:\n"
               "%d - %d is not within 0 - %d\n", x, x+w, width);
        return Iter(base + y*ystride + t*tstride + c*cstride, width);
    }
    #else
    struct Iter {
        const float *const addr;
        Iter() : addr(NULL) {}
        Iter(const float *a) : addr(a) {}
        float operator[](int x) const {return addr[x];}
        ImageStack::Expr::Vec::type vec(int x) const {
            return ImageStack::Expr::Vec::load(addr+x);
        }
    };
    Iter scanline(int x, int y, int t, int c, int w) const {
        return Iter(base + y*ystride + t*tstride + c*cstride);
    }
    #endif
    
    bool boundedVecX() const {return false;}
    int minVecX() const {return -Expr::HUGE_INT;}
    int maxVecX() const {return Expr::HUGE_INT;}

    void prepare(Expr::Region r, int phase) const {
        assert(r.x >= 0 && r.x+r.width <= width &&
               r.y >= 0 && r.y+r.height <= height &&
               r.t >= 0 && r.t+r.frames <= frames &&
               r.c >= 0 && r.c+r.channels <= channels, 
               "Expression would access image out of bounds: %d %d %d %d  %d %d %d %d\n",
               r.x, r.y, r.t, r.c, r.width, r.height, r.frames, r.channels);
    }

    std::pair<float, float> bounds(Expr::Region r) const {
        return make_pair(-INF, INF);
    }


    // Construct an image from a bounded expression object. 
    template<typename T>
    Image(const T &expr_, const FloatExprType(T) *ptr = NULL) :
        width(0), height(0), frames(0), channels(0),
        ystride(0), tstride(0), cstride(0), data(), base(NULL) {
        FloatExprType(T) expr(expr_);
        assert(expr.getSize(0) && expr.getSize(1) && expr.getSize(2) && expr.getSize(3),
               "Can only construct an image from a bounded expression\n");
        (*this) = Image(expr.getSize(0), expr.getSize(1), expr.getSize(2), expr.getSize(3));
        set(expr);
    }

    Image(const Image &other) :
        width(other.width), height(other.height), frames(other.frames), channels(other.channels),
        ystride(other.ystride), tstride(other.tstride), cstride(other.cstride),
        data(other.data), base(other.base) {
    }

private:


    template<int outChannels, typename A, typename B, typename C, typename D>
    void setChannelsGeneric(const A &exprA,
                            const B &exprB,
                            const C &exprC,
                            const D &exprD) const {
        int wA = exprA.getSize(0), hA = exprA.getSize(1), fA = exprA.getSize(2);
        int wB = exprB.getSize(0), hB = exprB.getSize(1), fB = exprB.getSize(2);
        int wC = exprC.getSize(0), hC = exprC.getSize(1), fC = exprC.getSize(2);
        int wD = exprD.getSize(0), hD = exprD.getSize(1), fD = exprD.getSize(2);
        assert(channels == outChannels,
               "The number of channels must equal the number of arguments\n");
        assert(exprA.getSize(3) <= 1 &&
               exprB.getSize(3) <= 1 &&
               exprC.getSize(3) <= 1 &&
               exprD.getSize(3) <= 1,
               "Each argument must be unbounded across channels or single-channel\n");
        assert((width == wA || wA == 0) &&
               (height == hA || hA == 0) &&
               (frames == fA || fA == 0),
               "Can only assign from sources of matching size\n");
        assert((width == wB || wB == 0) &&
               (height == hB || hB == 0) &&
               (frames == fB || fB == 0),
               "Can only assign from sources of matching size\n");
        assert((width == wC || wC == 0) &&
               (height == hC || hC == 0) &&
               (frames == fC || fC == 0),
               "Can only assign from sources of matching size\n");
        assert((width == wD || wD == 0) &&
               (height == hD || hD == 0) &&
               (frames == fD || fD == 0),
               "Can only assign from sources of matching size\n");

        bool boundedVX = (exprA.boundedVecX() || 
                          exprB.boundedVecX() || 
                          exprC.boundedVecX() || 
                          exprD.boundedVecX());       
        int minVX = max(max(exprA.minVecX(), exprB.minVecX()),
                        max(exprC.minVecX(), exprD.minVecX()));
        int maxVX = min(min(exprA.maxVecX(), exprB.maxVecX()),
                        min(exprC.maxVecX(), exprD.maxVecX()));

        Expr::Region r = {0, 0, 0, 0, width, height, frames, 1};
        
        exprA.prepare(r, 0);
        exprB.prepare(r, 0);
        exprC.prepare(r, 0);
        exprD.prepare(r, 0);
        exprA.prepare(r, 1);
        exprB.prepare(r, 1);
        exprC.prepare(r, 1);
        exprD.prepare(r, 1);
        exprA.prepare(r, 2);
        exprB.prepare(r, 2);
        exprC.prepare(r, 2);
        exprD.prepare(r, 2);

        // 4 or 8-wide vector code, distributed across cores
        for (int t = 0; t < frames; t++) {

            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (int y = 0; y < height; y++) {
                const int w = width;
                const int cs = cstride;

                const typename A::Iter iterA = exprA.scanline(0, y, t, 0, w);
                const typename B::Iter iterB = exprB.scanline(0, y, t, 0, w);
                const typename C::Iter iterC = exprC.scanline(0, y, t, 0, w);
                const typename D::Iter iterD = exprD.scanline(0, y, t, 0, w);

                float *const dst1 = base + t*tstride + y*ystride;
                float *const dst2 = outChannels > 1 ? dst1 + cs : NULL;
                float *const dst3 = outChannels > 2 ? dst2 + cs : NULL;
                float *const dst4 = outChannels > 3 ? dst3 + cs : NULL;

                Expr::setScanlineMulti(iterA, iterB, iterC, iterD,
                                       dst1, dst2, dst3, dst4,
                                       0, width, 
                                       boundedVX, minVX, maxVX);                
            }
        }

        exprA.prepare(r, 3);
        exprB.prepare(r, 3);
        exprC.prepare(r, 3);
        exprD.prepare(r, 3);
    }




    struct Payload {
        Payload(size_t size) : data(NULL) {
            // In some cases we don't need to clear the memory, but
            // typically this is optimized away by the system, so we
            // don't care. On linux it just mmaps /dev/zero.
            data = (float *)calloc(size, sizeof(float));
            //printf("Allocating %d bytes\n", (int)(size * sizeof(float)));
            if (!data) {
                panic("Could not allocate %d bytes for image data\n",
                      size * sizeof(float));
            }
        }
        ~Payload() {
            free(data);
        }
        float *data;
    private:
        // These are private to prevent copying a Payload
        Payload(const Payload &other) : data(NULL) {}
        void operator=(const Payload &other) {data = NULL;}
    };

    // Compute a 32-byte aligned address within data
    static float *compute_base(const std::shared_ptr<const Payload> &payload) {
        float *base = payload->data;
        while (((size_t)base) & 0x1f) base++;
        return base;
    }

    // Region constructor
    Image(const Image &im, int x, int y, int t, int c,
          int xs, int ys, int ts, int cs) :
        width(xs), height(ys), frames(ts), channels(cs),
        ystride(im.ystride), tstride(im.tstride), cstride(im.cstride),
        data(im.data), base(&im(x, y, t, c)) {
        // Note that base is no longer aligned. You're only guaranteed
        // alignment if you allocate an image from scratch.
        assert(xs > 0 && ys > 0 && ts > 0 && cs > 0,
               "Region must have strictly positive size: %d %d %d %d\n", xs, ys, ts, cs);
        /* Some regions may venture into invalid areas to help with indexing math 
        assert(x >= 0 && x+xs <= im.width &&
               y >= 0 && y+ys <= im.height &&
               t >= 0 && t+ts <= im.frames &&
               c >= 0 && c+cs <= im.channels,
               "Region must fit within original image\n");
        */
    }

    std::shared_ptr<const Payload> data;
    float *base;
};

template<typename SX, typename SY, typename ST, typename SC, bool AffineCase, bool ShiftedCase>
struct ImRefIter;


// Iterate across a scanline of a function that we're sampling in an unrestricted manner
template<typename SX, typename SY, typename ST, typename SC>
class ImRefIter<SX, SY, ST, SC, false, false> {
    Image im;
    const typename SX::Iter sx;
    const typename SY::Iter sy;
    const typename ST::Iter st;
    const typename SC::Iter sc;
public:
    ImRefIter() {}
    ImRefIter(Image im_,
              const typename SX::Iter &sx_,
              const typename SY::Iter &sy_,
              const typename ST::Iter &st_,
              const typename SC::Iter &sc_) : 
        im(im_), 
        sx(sx_), sy(sy_), st(st_), sc(sc_) {
    }
    float operator[](int x) const {
        return im(sx[x], sy[x], st[x], sc[x]);
    }
    Vec::type vec(int x) const {                
        union {
            float f[Vec::width];
            Vec::type v;
        } v;
        for (int i = 0; i < Vec::width; i++) {
            v.f[i] = (*this)[x+i];                   
        }
        return v.v;
    }                                            
};
    
// Iterate across a scanline of an image where
// 1) The index in X is affine
// 2) No other indices depend on X
template<typename SX, typename SY, typename ST, typename SC>
class ImRefIter<SX, SY, ST, SC, true, false> {
    Expr::AffineSampleX<Image>::Iter iter;
public:
    ImRefIter() {}
    ImRefIter(Image im, 
              const typename SX::Iter &sx,
              const typename SY::Iter &sy,
              const typename ST::Iter &st,
              const typename SC::Iter &sc) : 
        iter(im.scanline(0, sy[0], st[0], sc[0], im.width), 
             sx[1] - sx[0], sx[0]) {
    }
    float operator[](int x) const {
        return iter[x];
    }
    Vec::type vec(int x) const {
        return iter.vec(x);
    }                             
};
    
// Iterate across a scanline of an image where
// 1) The index in X is X + constant
// 2) No other indices depend on X
template<typename SX, typename SY, typename ST, typename SC>
class ImRefIter<SX, SY, ST, SC, true, true> {
    Expr::_Shift<Image>::Iter iter;
public:
    ImRefIter() {}
    ImRefIter(Image im, 
              const typename SX::Iter &sx,
              const typename SY::Iter &sy,
              const typename ST::Iter &st,
              const typename SC::Iter &sc) : 
        iter(im.scanline(0, sy[0], st[0], sc[0], im.width), -sx[0]) {
    }
    float operator[](int x) const {
        return iter[x];
    }
    Vec::type vec(int x) const {
        return iter.vec(x);
    }                             
};
    
// A computed reference to an image site (e.g. im(X+1, Y, sin(T), 0))
template<typename SX, typename SY, typename ST, typename SC, bool AffineCase, bool ShiftedCase>
class ImageRef {
    const Image im;
    const SX sx;
    const SY sy;
    const ST st;
    const SC sc;            
    int sizes[4];
public:
    typedef ImageRef<SX, SY, ST, SC, AffineCase, ShiftedCase> FloatExpr;
        
    static const bool dependsOnX = true;
        
    ImageRef(Image im_,
             const SX &sx_, 
             const SY &sy_, 
             const ST &st_, 
             const SC &sc_) : 
        im(im_), sx(sx_), sy(sy_), st(st_), sc(sc_) {
            
        for (int i = 0; i < 4; i++) {                    
            sizes[i] = std::max(std::max(sx.getSize(i), sy.getSize(i)), 
                                std::max(st.getSize(i), sc.getSize(i)));
            assert(sx.getSize(i) == 0 || sx.getSize(i) == sizes[i], 
                   "X coordinate must be unbounded or have the same size as other coordinates\n");
            assert(sy.getSize(i) == 0 || sy.getSize(i) == sizes[i], 
                   "Y coordinate must be unbounded or have the same size as other coordinates\n");
            assert(st.getSize(i) == 0 || st.getSize(i) == sizes[i], 
                   "T coordinate must be unbounded or have the same size as other coordinates\n");
            assert(sc.getSize(i) == 0 || sc.getSize(i) == sizes[i], 
                   "C coordinate must be unbounded or have the same size as other coordinates\n");
        }
    }

    int getSize(int i) const {
        return sizes[i];
    }

    typedef ImRefIter<SX, SY, ST, SC, AffineCase, ShiftedCase> Iter;
        
    Iter scanline(int x, int y, int t, int c, int width) const {
        return Iter(im, 
                    sx.scanline(x, y, t, c, width),
                    sy.scanline(x, y, t, c, width),
                    st.scanline(x, y, t, c, width),
                    sc.scanline(x, y, t, c, width));
    }

    // The image is safely over-allocated so that you can always pull vectors from it
    bool boundedVecX() const {
        return false;
    }
    int minVecX() const {
        return -Expr::HUGE_INT;
    }
    int maxVecX() const {
        return Expr::HUGE_INT;
    }

    std::pair<float, float> bounds(Expr::Region r) const {
        return make_pair(-INF, INF);
    }

    void prepare(Expr::Region r, int phase) const {
        // Figure out what the arguments require
        sx.prepare(r, phase);
        sy.prepare(r, phase);
        st.prepare(r, phase);
        sc.prepare(r, phase);

        // Prepare the image, which just amounts to a bounds check
        std::pair<int, int> xb = sx.bounds(r);
        std::pair<int, int> yb = sy.bounds(r);
        std::pair<int, int> tb = st.bounds(r);
        std::pair<int, int> cb = sc.bounds(r);
            
        Expr::Region r2 = {xb.first, yb.first, tb.first, cb.first, 
                           xb.second - xb.first + 1,
                           yb.second - yb.first + 1,
                           tb.second - tb.first + 1,
                           cb.second - cb.first + 1};
        im.prepare(r2, phase);
    }
};

// Sample an image at a computed address. Returns a different
// type based on static analysis of the arguments.
template<typename SX, typename SY, typename ST, typename SC>
ImageRef<IntExprType(SX), IntExprType(SY), IntExprType(ST), IntExprType(SC), 
         AffineCase(SX, SY, ST, SC), 
         ShiftedCase(SX, SY, ST, SC)>
Image::operator()(const SX &x, const SY &y, const ST &t, const SC &c) const {
    return ImageRef<IntExprType(SX), IntExprType(SY), IntExprType(ST), IntExprType(SC), 
                    AffineCase(SX, SY, ST, SC), 
                    ShiftedCase(SX, SY, ST, SC)>
    ((*this), x, y, t, c);
}

template<typename SX, typename SY, typename SC>
ImageRef<IntExprType(SX), IntExprType(SY), Expr::ConstInt, IntExprType(SC), 
         AffineCase(SX, SY, Expr::ConstInt, SC), ShiftedCase(SX, SY, Expr::ConstInt, SC)>
Image::operator()(const SX &x, const SY &y, const SC &c) const {
    return ImageRef<IntExprType(SX), IntExprType(SY), Expr::ConstInt, IntExprType(SC), 
                    AffineCase(SX, SY, Expr::ConstInt, SC), 
                    ShiftedCase(SX, SY, Expr::ConstInt, SC)>
    ((*this), x, y, 0, c);
}

template<typename SX, typename SY>
ImageRef<IntExprType(SX), IntExprType(SY), Expr::ConstInt, Expr::ConstInt, 
         AffineCase(SX, SY, Expr::ConstInt, Expr::ConstInt), 
         ShiftedCase(SX, SY, Expr::ConstInt, Expr::ConstInt)>
Image::operator()(const SX &x, const SY &y) const {
    return ImageRef<IntExprType(SX), IntExprType(SY), Expr::ConstInt, Expr::ConstInt, 
                    AffineCase(SX, SY, Expr::ConstInt, Expr::ConstInt), 
                    ShiftedCase(SX, SY, Expr::ConstInt, Expr::ConstInt)>
    ((*this), x, y, 0, 0);
}

template<typename SX>
ImageRef<IntExprType(SX), Expr::ConstInt, Expr::ConstInt, Expr::ConstInt, 
         AffineCase(SX, Expr::ConstInt, Expr::ConstInt, Expr::ConstInt), 
         ShiftedCase(SX, Expr::ConstInt, Expr::ConstInt, Expr::ConstInt)>
Image::operator()(const SX &x) const {
    return ImageRef<IntExprType(SX), Expr::ConstInt, Expr::ConstInt, Expr::ConstInt, 
                    AffineCase(SX, Expr::ConstInt, Expr::ConstInt, Expr::ConstInt), 
                    ShiftedCase(SX, Expr::ConstInt, Expr::ConstInt, Expr::ConstInt)>
    ((*this), x, 0, 0, 0);
}

#undef AffineCase        
#undef ShiftedCase


#include "footer.h"
#endif
