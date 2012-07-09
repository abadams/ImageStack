#ifndef IMAGESTACK_IMAGE_H
#define IMAGESTACK_IMAGE_H

#include "Lazy.h"

#include "tables.h"
#include "header.h"

// The image data type.

// It's a reference-counted pointer type.

// Note that "const Image" means that the reference doesn't change, not
// that the pixel data doesn't. It's equivalent to "float * const
// foo", not "const float * foo". This means that methods tagged const
// are those which do not change the metadata, not those which do not
// change the pixel data.

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
        ystride(w), tstride(w *h), cstride(w *h *f),
        data(new Payload(w *h *f *c+7)), base(compute_base(data)) {
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
        return region(0, x, 0, 0, width, x, frames, channels);
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

    void operator+=(const float f) const {
        set((*this) + f);
    }

    void operator*=(const float f) const {
        set((*this) * f);
    }

    void operator-=(const float f) const {
        set((*this) - f);
    }

    void operator/=(const float f) const {
        set((*this) / f);
    }

    template<typename A, typename B, typename Enable = typename A::Lazy>
    struct LazyCheck {
        typedef B result;
    };

    template<typename T>
    typename LazyCheck<T, void>::result operator+=(const T &other) const {
        set((*this) + other);
    }

    template<typename T>
    typename LazyCheck<T, void>::result operator*=(const T &other) const {
        set((*this) * other);
    }

    template<typename T>
    typename LazyCheck<T, void>::result operator-=(const T &other) const {
        set((*this) - other);
    }

    template<typename T>
    typename LazyCheck<T, void>::result operator/=(const T &other) const {
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

    // A macro to check if a type is castable to a lazy expression type
    #define LazyType(T) typename ImageStack::Lazy::Lazyable<T>::t

    // Evaluate a expression object defined in Lazy.h
    // The second template argument prevents instantiations from
    // things that don't satisfy the trait "lazyable"
    template<typename T>
    void set(const T expr_, const LazyType(T) *enable = NULL) const {
        LazyType(T) expr(expr_);
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

        for (int c = 0; c < channels; c++) {
            for (int t = 0; t < frames; t++) {
                for (int yt = 0; yt < height; yt += 16) {
                    const int maxY = min(yt + 16, height);
                    for (int xt = 0; xt < width; xt += 1024) {
                        const int maxX = min(xt + 1024, width);
                        // Assume we chunk here for now. This decision assumes
                        // no footprint in C and T, and could be grossly
                        // wasteful! In future, we should just explicitly
                        // check the footprint and use that to decide where to
                        // chunk.  In fact, with the current logic inside
                        // Func, we only keep track of which scanlines were
                        // evaluated, not which channel or frames, so we're
                        // going to produce bogus output. 
                        expr.prepare(0, xt, yt, t, c, maxX-xt, maxY-yt, 1, 1);
                        expr.prepare(1, xt, yt, t, c, maxX-xt, maxY-yt, 1, 1);

                        #ifdef _OPENMP
                        #pragma omp parallel for
                        #endif
                        for (int y = yt; y < maxY; y++) {
                            LazyType(T)::Iter iter = expr.scanline(xt, y, t, c, maxX-xt);
                            float *const dst = base + c*cstride + t*tstride + y*ystride;
                            ImageStack::Lazy::setScanline(iter, dst, xt, maxX, boundedVX, minVX, maxVX);
                        }                        
                    }
                }
            }
        }
        // Clean up any resources
        expr.prepare(2, 0, 0, 0, 0, width, height, frames, channels);
    }


    // A version of set that takes a set of up to 4 expressions and
    // sets the image's channels accordingly. This is more efficient
    // than calling set repeatedly when the expressions share common
    // subexpressions. At each pixel, each expression is evaluated,
    // and then each value is stored, preventing nasty chicken-and-egg
    // problems when, for example, permuting channels.
    template<typename A, typename B, typename C, typename D>
    void setChannels(const A a, const B b, const C c, const D d,
                     LazyType(A) *pa = NULL,
                     LazyType(B) *pb = NULL,
                     LazyType(C) *pc = NULL,
                     LazyType(D) *pd = NULL) const {
        setChannelsGeneric<4, LazyType(A), LazyType(B), LazyType(C), LazyType(D)>(
            LazyType(A)(a),
            LazyType(B)(b),
            LazyType(C)(c),
            LazyType(D)(d));
    }

    template<typename A, typename B, typename C>
    void setChannels(const A &a, const B &b, const C &c,
                     LazyType(A) *pa = NULL,
                     LazyType(B) *pb = NULL,
                     LazyType(C) *pc = NULL) const {
        setChannelsGeneric<3, LazyType(A), LazyType(B), LazyType(C), LazyType(float)>(
            LazyType(A)(a),
            LazyType(B)(b),
            LazyType(C)(c),
            LazyType(float)(0));
    }

    template<typename A, typename B>
    void setChannels(const A &a, const B &b,
                     LazyType(A) *pa = NULL,
                     LazyType(B) *pb = NULL) const {
        setChannelsGeneric<2, LazyType(A), LazyType(B), LazyType(float), LazyType(float)>(
            LazyType(A)(a),
            LazyType(B)(b),
            LazyType(float)(0),
            LazyType(float)(0));
    }

    // An image itself is one such expression thing. Here's the
    // interface it needs to implement for that to happen:
    typedef Image Lazy;
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
        ImageStack::Lazy::Vec::type vec(int x) const {
            assert(x >= 0 && x <= width - ImageStack::Lazy::Vec::width,
                   "Vector access out of bounds in image iterator:\n"
                   "%d is not sufficiently within 0 - %d\n", x, width);
            return ImageStack::Lazy::Vec::load(addr+x);
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
        ImageStack::Lazy::Vec::type vec(int x) const {
            return ImageStack::Lazy::Vec::load(addr+x);
        }
    };
    Iter scanline(int x, int y, int t, int c, int w) const {
        return Iter(base + y*ystride + t*tstride + c*cstride);
    }
    #endif

    bool boundedVecX() const {return true;}
    int minVecX() const {return 0;}
    int maxVecX() const {return width-ImageStack::Lazy::Vec::width;}

    void prepare(int phase, int x, int y, int t, int c,
                 int xs, int ys, int ts, int cs) const {
        assert(x >= 0 && x+xs <= width &&
               y >= 0 && y+ys <= height &&
               t >= 0 && t+ts <= frames &&
               c >= 0 && c+cs <= channels, 
               "Expression would access image out of bounds: %d %d %d %d  %d %d %d %d\n",
               x, y, t, c, xs, ys, ts, cs);
    }

    // Construct an image from a bounded expression object. No consts
    // allowed, so we require a ::Lazy subtype instead of the more
    // general LazyType(T) macro.
    template<typename T>
    Image(const T &func, const typename T::Lazy *ptr = NULL) :
        width(0), height(0), frames(0), channels(0),
        ystride(0), tstride(0), cstride(0), data(), base(NULL) {
        assert(func.getSize(0) && func.getSize(1) && func.getSize(2) && func.getSize(3),
               "Can only construct an image from a bounded expression\n");
        (*this) = Image(func.getSize(0), func.getSize(1), func.getSize(2), func.getSize(3));
        set(func);
    }

    Image(const Image &other) :
        width(other.width), height(other.height), frames(other.frames), channels(other.channels),
        ystride(other.ystride), tstride(other.tstride), cstride(other.cstride),
        data(other.data), base(other.base) {
    }

private:


    template<int outChannels, typename A, typename B, typename C, typename D>
    void setChannelsGeneric(const A &funcA,
                            const B &funcB,
                            const C &funcC,
                            const D &funcD) const {
        int wA = funcA.getSize(0), hA = funcA.getSize(1), fA = funcA.getSize(2);
        int wB = funcB.getSize(0), hB = funcB.getSize(1), fB = funcB.getSize(2);
        int wC = funcC.getSize(0), hC = funcC.getSize(1), fC = funcC.getSize(2);
        int wD = funcD.getSize(0), hD = funcD.getSize(1), fD = funcD.getSize(2);
        assert(channels == outChannels,
               "The number of channels must equal the number of arguments\n");
        assert(funcA.getSize(3) <= 1 &&
               funcB.getSize(3) <= 1 &&
               funcC.getSize(3) <= 1 &&
               funcD.getSize(3) <= 1,
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

        // TODO: respect minVecX and maxVecX

        // 4 or 8-wide vector code, distributed across cores
        const int vec_width = ImageStack::Lazy::Vec::width;
        for (int t = 0; t < frames; t++) {


            funcA.prepare(0, 0, 0, t, 0, width, height, 1, 1);
            funcB.prepare(0, 0, 0, t, 0, width, height, 1, 1);
            funcC.prepare(0, 0, 0, t, 0, width, height, 1, 1);
            funcD.prepare(0, 0, 0, t, 0, width, height, 1, 1);
            funcA.prepare(1, 0, 0, t, 0, width, height, 1, 1);
            funcB.prepare(1, 0, 0, t, 0, width, height, 1, 1);
            funcC.prepare(1, 0, 0, t, 0, width, height, 1, 1);
            funcD.prepare(1, 0, 0, t, 0, width, height, 1, 1);
            #ifdef _OPENMP
            #pragma omp parallel for
            #endif
            for (int y = 0; y < height; y++) {
                const int w = width;
                const int cs = cstride;
                float *const dst = base + t*tstride + y*ystride;
                const typename A::Iter iterA = funcA.scanline(0, y, t, 0, w);
                const typename B::Iter iterB = funcB.scanline(0, y, t, 0, w);
                const typename C::Iter iterC = funcC.scanline(0, y, t, 0, w);
                const typename D::Iter iterD = funcD.scanline(0, y, t, 0, w);

                int x = 0;

                if (vec_width > 1 && w > vec_width*2) {
                    while (x < (w-(ImageStack::Lazy::Vec::width-1))) {
                        ImageStack::Lazy::Vec::type va, vb, vc, vd;

                        // Compute the value of each input at this pixel
                        va = iterA.vec(x);
                        if (outChannels > 1) vb = iterB.vec(x);
                        if (outChannels > 2) vc = iterC.vec(x);
                        if (outChannels > 3) vd = iterD.vec(x);

                        // Store the results
                        ImageStack::Lazy::Vec::store(va, dst + x);
                        if (outChannels > 1)
                            ImageStack::Lazy::Vec::store(vb, dst + cs + x);
                        if (outChannels > 2)
                            ImageStack::Lazy::Vec::store(vc, dst + cs*2 + x);
                        if (outChannels > 3)
                            ImageStack::Lazy::Vec::store(vd, dst + cs*3 + x);
                        x += vec_width;
                    }
                }

                // Scalar part at the end
                while (x < w) {
                    float va, vb, vc, vd;

                    va = iterA[x];
                    if (outChannels > 1) vb = iterB[x];
                    if (outChannels > 2) vc = iterC[x];
                    if (outChannels > 3) vd = iterD[x];

                    dst[x] = va;
                    if (outChannels > 1)
                        dst[x + cs] = vb;
                    if (outChannels > 2)
                        dst[x + cs*2] = vc;
                    if (outChannels > 3)
                        dst[x + cs*3] = vd;
                    x++;
                }
            }
        }
        funcA.prepare(2, 0, 0, 0, 0, width, height, frames, 1);
        funcB.prepare(2, 0, 0, 0, 0, width, height, frames, 1);
        funcC.prepare(2, 0, 0, 0, 0, width, height, frames, 1);
        funcD.prepare(2, 0, 0, 0, 0, width, height, frames, 1);
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

// Clean up after myself
#undef LazyType




#include "footer.h"
#endif
