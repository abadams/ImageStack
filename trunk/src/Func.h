#ifndef IMAGESTACK_FUNC_H
#define IMAGESTACK_FUNC_H

#include "Image.h"
#include "header.h"

namespace Expr {

    // Is an expression of the form a x + b
    #define AffineInX(T) _AffineInX<T>::value

    template<typename T>
    struct _AffineInX {
        static const bool value = false;
    };
    
    template<>
    struct _AffineInX<X> {
        static const bool value = true;
    };
    
    template<>
    struct _AffineInX<int> {
        static const bool value = true;
    };

    template<>
    struct _AffineInX<ConstInt> {
        static const bool value = true;
    };
    
    template<typename A, typename B>
    struct _AffineInX<IBinaryOp<A, B, Vec::Add> > {
        static const bool value = AffineInX(A) && AffineInX(B);
    };
    
    template<typename A, typename B>
    struct _AffineInX<IBinaryOp<A, B, Vec::Sub> > {
        static const bool value = AffineInX(A) && AffineInX(B);
    };
    
    template<typename A>
    struct _AffineInX<IBinaryOp<ConstInt, A, Vec::Mul> > {
        static const bool value = AffineInX(A);
    };
    
    template<typename A>
    struct _AffineInX<IBinaryOp<A, ConstInt, Vec::Mul> > {
        static const bool value = AffineInX(A);
    };
    
    // Is an expression of the form x + b
    #define ShiftedInX(T) _ShiftedInX<T>::value

    template<typename T>
    struct _ShiftedInX {
        static const bool value = false;
    };

    template<>
    struct _ShiftedInX<X> {
        static const bool value = true;
    };

    template<typename A>
    struct _ShiftedInX<IBinaryOp<A, ConstInt, Vec::Add> > {
        static const bool value = ShiftedInX(A);
    };

    template<typename A>
    struct _ShiftedInX<IBinaryOp<A, ConstInt, Vec::Sub> > {
        static const bool value = ShiftedInX(A);
    };

    // Is an expression constant per scanline
    #define IndependentOfX(T) _IndependentOfX<T>::value

    template<typename T>
    struct _IndependentOfX {
        static const bool value = !T::dependsOnX;
    };

    template<>
    struct _IndependentOfX<int> {
        static const bool value = true;
    };

    template<>
    struct _IndependentOfX<float> {
        static const bool value = true;
    };

    struct BaseFunc {
        Image im;
        int minX, minY, minT, minC;
        int maxX, maxY, maxT, maxC;
        bool lazy;

        vector<int> evaluated;

        // Evaluate myself at the given scanline only if necessary
        void evalScanlineIfNeeded(int y, int t, int c) {
            if (!lazy) return;
           
            int idx = ((c-minC) * im.frames + t-minT) * im.height + y-minY;
            if (!evaluated[idx])  {
                // TODO: consider locking the scanline during
                // evaluation. As it stands, if multiple threads
                // hammer on the same scanline, there will be wasted
                // work (and probably lots of fighting over cache
                // lines).  Fortunately, that's not how openmp
                // schedules its threads.
                evalScanline(y, t, c);
                evaluated[idx] = true;
            }
            
        }

        _Shift<Image>::Iter scanline(int x, int y, int t, int c, int width) {
            evalScanlineIfNeeded(y, t, c);
            Image::Iter iter = im.scanline(x-minX, y-minY, t-minT, c-minC, width);
            return _Shift<Image>::Iter(iter, minX);            
        }

        // Evaluate myself into buffer at the given scanline. 
        virtual void evalScanline(int y, int t, int c) = 0;

        // Prepare to be evaluated over a given region
        virtual void prepareFunc(int phase, Region r) = 0;

        // Evaluate your expression into an image (i.e. call im.set(expr))
        virtual void realize(Image im) = 0;

        virtual int getSize(int i) = 0;
    };

    template<typename T, typename Enable = typename T::FloatExpr>
    struct DerivedFunc : public BaseFunc {
        const T expr;
        int lastPhase;
        int minVecX;
        int maxVecX;
        bool boundedVecX;

        DerivedFunc(const T &e) : 
            expr(e), lastPhase(-1), 
            minVecX(e.minVecX()), maxVecX(e.maxVecX()), boundedVecX(e.boundedVecX()) {}

        void realize(Image m) {
            m.set(expr);
        }

        void evalScanline(int y, int t, int c) {
            //printf("Evaluating %p at scanline %d-%d %d %d %d\n",
            //this, minX, maxX, y, t, c);
            //printf("Image has size: %d %d %d %d\n", im.width, im.height, im.frames, im.channels);
            //printf("Computing destination address...\n");
            float *const dst = &im(0, y-minY, t-minT, c-minC) - minX;

            //printf("Destination address is %p.\n"
            //"Computing source iterator...\n", dst);
            typename T::Iter src = expr.scanline(minX, y, t, c, maxX-minX);
            
            //printf("Done. Running kernel...\n");

            setScanline(src, dst, minX, maxX, 
                        boundedVecX, minVecX, maxVecX);

            //printf("Done\n"); fflush(stdout);
        }   

        void prepareFunc(int phase, Region r) {

            // recurse
            expr.prepare(phase, r);

            //printf("Preparing %p at phase %d over region %d %d %d %d  %d %d %d %d\n", this,
            //phase, r.x, r.y, r.t, r.c, r.width, r.height, r.frames, r.channels);               

            // Each phase we get called multiple times according to
            // how many times we occur in the expression

            // The start of a new page
            // In phase zero we take unions of the regions requested
            if (phase == 0) {
                if (lastPhase != 0) {
                    minX = r.x;
                    minY = r.y;
                    minT = r.t;
                    minC = r.c;
                    maxX = r.x + r.width;
                    maxY = r.y + r.height;
                    maxT = r.t + r.frames;
                    maxC = r.c + r.channels;
                } else {
                    minX = std::min(minX, r.x);
                    minY = std::min(minY, r.y);
                    minT = std::min(minT, r.t);
                    minC = std::min(minC, r.c);
                    maxX = std::max(maxX, r.x + r.width);
                    maxY = std::max(maxY, r.y + r.height);
                    maxT = std::max(maxT, r.t + r.frames);
                    maxC = std::max(maxC, r.c + r.channels);
                }
                //printf("New bounds are %d %d %d %d - %d %d %d %d\n", 
                //minX, minY, minT, minC, maxX-minX, maxY-minY, maxT-minT, maxC-minC);
            } else if (phase == 1) {
                // Still at the start of a new page
                // In phase one (the first time we are called), we prepare storage
                if (lastPhase != 1) {
                    if (!im.defined() ||                            
                        im.width < maxX - minX ||
                        im.height < maxY - minY ||
                        im.frames < maxT - minT ||
                        im.channels < maxC - minC) {
                        //printf("Allocating backing for %p %d %d %d %d\n", 
                        //this, maxX-minX, maxY-minY, maxT-minT, maxC-minC);
                        im = Image(maxX - minX, maxY - minY, maxT - minT, maxC - minC);
                        //printf("Done\n");

                    } else {
                        // no need to allocate
                    }

                    if (lazy) {
                        // No scanlines have been evaluated
                        evaluated.assign((maxY - minY)*(maxT - minT)*(maxC - minC), false);
                    } else {
                        // evaluate all scanlines here
                        for (int ec = minC; ec < maxC; ec++) {
                            for (int et = minT; et < maxT; et++) {
                                #ifdef _OPENMP
                                #pragma omp parallel for
                                #endif
                                for (int ey = minY; ey < maxY; ey++) {
                                    evalScanline(ey, et, ec);
                                }
                            }
                        }
                    }                    
                }
            } else if (phase == 2) {
                if (lastPhase != 2) {
                    // Clean-up
                    // We're done iterating over everything. Clean up.
                    im = Image();
                    //printf("Freeing backing for %p\n", this);
                }
            }       

            lastPhase = phase;

        }

        int getSize(int i) {
            return expr.getSize(i);
        }
    };

    class Func {
        std::shared_ptr<BaseFunc> ptr;
    public:

        template<typename T>
        Func(const T t, FloatExprType(T) *enable = NULL) {            
            ptr.reset(new DerivedFunc<FloatExprType(T)>(t));
            ptr->lazy = true;
        }

        Func() {
        }

        // Implements the float expr interface
        typedef Func FloatExpr;

        const static bool dependsOnX = true;

        template<typename SX, typename SY, typename ST, typename SC, bool AffineCase, bool ShiftedCase>
        struct FuncRefIter;
            
        // Iterate across a scanline of a function that we're sampling in an unrestricted manner
        template<typename SX, typename SY, typename ST, typename SC>
        class FuncRefIter<SX, SY, ST, SC, false, false> {
            const Image im;
            const int minX, minY, minT, minC;
            const typename SX::Iter sx;
            const typename SY::Iter sy;
            const typename ST::Iter st;
            const typename SC::Iter sc;
        public:
            FuncRefIter() {}
            FuncRefIter(const std::shared_ptr<BaseFunc> &f,
                        const typename SX::Iter &sx_,
                        const typename SY::Iter &sy_,
                        const typename ST::Iter &st_,
                        const typename SC::Iter &sc_) : 
                im(f->im), minX(f->minX), minY(f->minY), minT(f->minT), minC(f->minC),
                sx(sx_), sy(sy_), st(st_), sc(sc_) {
            }
            float operator[](int x) const {
                return im(sx[x]-minX, sy[x]-minY, st[x]-minT, sc[x]-minC);
            }
            Vec::type vec(int x) const {
                if (Vec::width == 8) {
                    return Vec::set((*this)[x],
                                    (*this)[x+1],
                                    (*this)[x+2],
                                    (*this)[x+3],
                                    (*this)[x+4],
                                    (*this)[x+5],
                                    (*this)[x+6],
                                    (*this)[x+7]);
                } else if (Vec::width == 4) {
                    return Vec::set((*this)[x],
                                    (*this)[x+1],
                                    (*this)[x+2],
                                    (*this)[x+3]);
                } else {
                    union {
                        float f[Vec::width];
                        Vec::type v;
                    } v;
                    for (int i = 0; i < Vec::width; i++) {
                        v.f[i] = (*this)[x];
                    }
                    return v.v;
                }
            }                                            
        };

        // Iterate across a scanline of a function where
        // 1) The index in X is affine
        // 2) No other indices depend on X
        template<typename SX, typename SY, typename ST, typename SC>
        class FuncRefIter<SX, SY, ST, SC, true, false> {
            AffineSampleX<Image>::Iter iter;
        public:
            FuncRefIter() {}
            FuncRefIter(const std::shared_ptr<BaseFunc> &f,
                        const typename SX::Iter &sx,
                        const typename SY::Iter &sy,
                        const typename ST::Iter &st,
                        const typename SC::Iter &sc) : 
                iter(f->im.scanline(0, sy[0]-f->minY, st[0]-f->minT, sc[0]-f->minC, f->im.width), 
                     sx[1] - sx[0], sx[0] - f->minX) {
                // Make sure the relevant scanline has been evaluated

                f->evalScanlineIfNeeded(sy[0], st[0], sc[0]);
            }
            float operator[](int x) const {
                return iter[x];
            }
            Vec::type vec(int x) const {
                return iter.vec(x);
            }                             
        };

        // Iterate across a scanline of a function where
        // 1) The index in X is X + constant
        // 2) No other indices depend on X
        template<typename SX, typename SY, typename ST, typename SC>
        class FuncRefIter<SX, SY, ST, SC, true, true> {
            _Shift<Image>::Iter iter;
        public:
            FuncRefIter() {}
            FuncRefIter(const std::shared_ptr<BaseFunc> &f,
                        const typename SX::Iter &sx,
                        const typename SY::Iter &sy,
                        const typename ST::Iter &st,
                        const typename SC::Iter &sc) : 
                iter(f->im.scanline(0, sy[0]-f->minY, st[0]-f->minT, sc[0]-f->minC, f->im.width), 
                     f->minX - sx[0]) {
                // Make sure the relevant scanline has been evaluated
                f->evalScanlineIfNeeded(sy[0], st[0], sc[0]);
            }
            float operator[](int x) const {
                return iter[x];
            }
            Vec::type vec(int x) const {
                return iter.vec(x);
            }                             
        };

        template<typename SX, typename SY, typename ST, typename SC, bool AffineCase, bool ShiftedCase>
        class FuncRef {
            std::shared_ptr<BaseFunc> ptr;
            const SX sx;
            const SY sy;
            const ST st;
            const SC sc;            
            int sizes[4];
        public:
            typedef FuncRef<SX, SY, ST, SC, AffineCase, ShiftedCase> FloatExpr;

            static const bool dependsOnX = true;
            
            FuncRef(const std::shared_ptr<BaseFunc> &ptr_,
                    const SX &sx_, 
                    const SY &sy_, 
                    const ST &st_, 
                    const SC &sc_) : 
                ptr(ptr_), sx(sx_), sy(sy_), st(st_), sc(sc_) {

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

                // If we're going to be sampling this func somewhat arbitrarily, then it can't be lazy
                if (!AffineCase && !ShiftedCase) ptr->lazy = false;
            }

            int getSize(int i) const {
                return sizes[i];
            }

            typedef FuncRefIter<SX, SY, ST, SC, AffineCase, ShiftedCase> Iter;

            Iter scanline(int x, int y, int t, int c, int width) const {
                return Iter(ptr,
                            sx.scanline(x, y, t, c, width),
                            sy.scanline(x, y, t, c, width),
                            st.scanline(x, y, t, c, width),
                            sc.scanline(x, y, t, c, width));
            }

            // In all cases we resolve to a load from an image, so
            // we're not particularly bounded in how we can vectorize
            bool boundedVecX() const {
                return false;
            }
            int minVecX() const {
                return 0xa0000000;
            }
            int maxVecX() const {
                return 0x3fffffff;
            }

            void prepare(int phase, Region r) const {
                // prepare the args over this range
                sx.prepare(phase, r);
                sy.prepare(phase, r);
                st.prepare(phase, r);
                sc.prepare(phase, r);

                // Prepare the func for more arbitrary sampling
                std::pair<int, int> xb = sx.bounds(r);
                std::pair<int, int> yb = sy.bounds(r);
                std::pair<int, int> tb = st.bounds(r);
                std::pair<int, int> cb = sc.bounds(r);

                // Make sure the region is adequately bounded
                Region r2 = {xb.first, yb.first, tb.first, cb.first, 
                             xb.second - xb.first + 1,
                             yb.second - yb.first + 1,
                             tb.second - tb.first + 1,
                             cb.second - cb.first + 1};
                ptr->prepareFunc(phase, r2);
            }
        };

#define ShiftedCase(SX, SY, ST, SC)                                     \
        (ShiftedInX(SX) && IndependentOfX(SY) && IndependentOfX(ST) && IndependentOfX(SC)) 
        
#define AffineCase(SX, SY, ST, SC)                                      \
        (AffineInX(SX) && IndependentOfX(SY) && IndependentOfX(ST) && IndependentOfX(SC)) 
        
        // Sample a function        
        template<typename SX, typename SY, typename ST, typename SC>
        FuncRef<IntExprType(SX), IntExprType(SY), IntExprType(ST), IntExprType(SC), 
                AffineCase(SX, SY, ST, SC), ShiftedCase(SX, SY, ST, SC)>
        operator()(const SX &x, const SY &y, const ST &t, const SC &c) const {
            return FuncRef<IntExprType(SX), IntExprType(SY), IntExprType(ST), IntExprType(SC), 
                           AffineCase(SX, SY, ST, SC), ShiftedCase(SX, SY, ST, SC)>
                (ptr, x, y, t, c);
        }

        template<typename SX, typename SY, typename SC>
        FuncRef<IntExprType(SX), IntExprType(SY), ConstInt, IntExprType(SC), 
                AffineCase(SX, SY, ConstInt, SC), ShiftedCase(SX, SY, ConstInt, SC)>
        operator()(const SX &x, const SY &y, const SC &c) const {
            return FuncRef<IntExprType(SX), IntExprType(SY), ConstInt, IntExprType(SC), 
                           AffineCase(SX, SY, ConstInt, SC), ShiftedCase(SX, SY, ConstInt, SC)>
                (ptr, x, y, 0, c);
        }

        template<typename SX, typename SY>
        FuncRef<IntExprType(SX), IntExprType(SY), ConstInt, ConstInt, 
                AffineCase(SX, SY, ConstInt, ConstInt), ShiftedCase(SX, SY, ConstInt, ConstInt)>
        operator()(const SX &x, const SY &y) const {
            return FuncRef<IntExprType(SX), IntExprType(SY), ConstInt, ConstInt, 
                           AffineCase(SX, SY, ConstInt, ConstInt), ShiftedCase(SX, SY, ConstInt, ConstInt)>
                (ptr, x, y, 0, 0);
        }

#undef AffineCase        
#undef ShiftedCase

        typedef _Shift<Image>::Iter Iter;
        Iter scanline(int x, int y, int t, int c, int width) const {        
            return ptr->scanline(x, y, t, c, width);
        }

        int getSize(int i) const {
            return ptr->getSize(i);
        }

        bool boundedVecX() const {return false;}
        int minVecX() const {return 0xa0000000;} 
        int maxVecX() const {return 0x3fffffff;}

        void prepare(int phase, Region r) const {
            ptr->prepareFunc(phase, r);
        }

        // Evaluate yourself into an existing image
        void realize(Image im) const {
            ptr->realize(im);
        }

        // Evaluate yourself into a new image
        Image realize(int w, int h, int f, int c) const {
            Image im(w, h, f, c);
            ptr->realize(im);
            return im;
        }

        // For funcs with bounds, evaluate over the bounds into a new
        // image
        Image realize() const {
            assert(getSize(0) && getSize(1) && getSize(2) && getSize(3), 
                   "Can only construct an image from a bounded function.\n"
                   "Call realize(width, height, frames, channels) instead.\n");
            Image im(getSize(0), getSize(1), getSize(2), getSize(3));
            ptr->realize(im);
            return im;
        }
    };

    inline void realizeFuncIntoImage(const Image &im, const Func &func) {    
        func.realize(im);
    }
    
    inline Image realizeFuncIntoNewImage(const Func &func) {
        return func.realize();
    }
}




#include "footer.h"
#endif
