#ifndef IMAGESTACK_FUNC_H
#define IMAGESTACK_FUNC_H

#include "Image.h"
#include "header.h"

namespace Expr {

    struct BaseFunc {
        Image im;
        int minX, minY, minT, minC;
        int maxX, maxY, maxT, maxC;
        bool lazy;

        vector<int> evaluated;

        // Evaluate myself into buffer at the given scanline. 
        virtual void evalScanline(int y, int t, int c) = 0;

        // Prepare to be evaluated over a given region
        virtual void prepareFunc(int phase, Region r) = 0;

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
                        // this, maxX-minX, maxY-minY, maxT-minT, maxC-minC);
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

    struct Func {
        std::shared_ptr<BaseFunc> ptr;

        template<typename T>
        Func(const T t, FloatExprType(T) *enable = NULL) {            
            ptr.reset(new DerivedFunc<FloatExprType(T)>(t));
            ptr->lazy = true;
        }

        Func() {
        }

        // Implements the float expr interface
        typedef Func FloatExpr;

        int getSize(int i) const {
            return ptr->getSize(i);
        }

        template<typename SX, typename SY, typename ST, typename SC>
        struct FuncRef {
            typedef FuncRef<SX, SY, ST, SC> FloatExpr;
            const SX sx;
            const SY sy;
            const ST st;
            const SC sc;
            std::shared_ptr<BaseFunc> ptr;
            FuncRef(const std::shared_ptr<BaseFunc> &ptr_,
                    const SX &sx_, 
                    const SY &sy_, 
                    const ST &st_, 
                    const SC &sc_) : 
                ptr(ptr_), sx(sx_), sy(sy_), st(st_), sc(sc_) {
                // We're going to be sampling this func somewhat arbitrarily, so it can't be lazy
                ptr->lazy = false;
            }

            int getSize(int i) {
                // TODO: max of sizes of sx, sy, st, sc
                return 0;
            }

            struct Iter {                
                const Image im;
                const int minX, minY, minT, minC;
                const typename SX::Iter sx;
                const typename SY::Iter sy;
                const typename ST::Iter st;
                const typename SC::Iter sc;
                Iter() {}
                Iter(const std::shared_ptr<BaseFunc> &f,
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
            Iter scanline(int x, int y, int t, int c, int width) const {
                return Iter(ptr,
                            sx.scanline(x, y, t, c, width),
                            sy.scanline(x, y, t, c, width),
                            st.scanline(x, y, t, c, width),
                            sc.scanline(x, y, t, c, width));                            
            }

            // We never call .vec on any dependents, so feel free to
            // vectorize this over any range - it resolves to a gather
            // from the backing image.
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

        template<typename SX, typename SY, typename ST, typename SC>
        FuncRef<IntExprType(SX), IntExprType(SY), IntExprType(ST), IntExprType(SC)>
        operator()(const SX &x, const SY &y, const ST &t, const SC &c) const {
            return FuncRef<IntExprType(SX), IntExprType(SY), IntExprType(ST), IntExprType(SC)>(ptr, x, y, t, c);
        }

        // TODO: special-case func samplings that return more restricted types

        typedef _Shift<Image>::Iter Iter;
        Iter scanline(int x, int y, int t, int c, int width) const {        
            int px = x-ptr->minX, py = y-ptr->minY, pt = t-ptr->minT, pc = c-ptr->minC;
            int idx = (pc * ptr->im.frames + pt) * ptr->im.height + py;
            if (ptr->lazy && !ptr->evaluated[idx])  {
                // TODO: consider locking the scanline during
                // evaluation. As it stands, if multiple threads
                // hammer on the same scanline, there will be wasted
                // work (and probably lots of fighting over cache
                // lines).  Fortunately, that's not how openmp
                // schedules its threads.
                ptr->evalScanline(y, t, c);
                ptr->evaluated[idx] = true;
            }
            Image::Iter iter = ptr->im.scanline(px, py, pt, pc, width);
            return _Shift<Image>::Iter(iter, ptr->minX);
        }

        bool boundedVecX() const {return false;}
        int minVecX() const {return 0xa0000000;} 
        int maxVecX() const {return 0x3fffffff;}

        void prepare(int phase, Region r) const {
            ptr->prepareFunc(phase, r);
        }
    };
}

#include "footer.h"
#endif
