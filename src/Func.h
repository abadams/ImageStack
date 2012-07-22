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
        std::string name;
        int size[4];

        vector<int> evaluated;

        // Evaluate myself at the given scanline only if necessary
        void evalScanlineIfNeeded(int y, int t, int c) {
            if (!lazy) return;
           
            int idx = ((c-minC) * im.frames + t-minT) * im.height + y;
            if (!evaluated[idx])  {
                // TODO: consider adding mem fences
                evaluated[idx] = 1;
                evalScanline(y, t, c);
                evaluated[idx] = 2;
            } else {            
                while (evaluated[idx] < 2) {
                    //printf("Stalling...\n");
                    sched_yield();
                }
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
        virtual void prepare(Region r, int phase) = 0;

        // Evaluate your expression into an image (i.e. call im.set(expr))
        virtual void realize(Image im) = 0;

        virtual std::pair<float, float> bounds(Region r) const = 0;

    };

    template<typename T, typename Enable = typename T::FloatExpr>
    struct DerivedFunc : public BaseFunc {
        const T expr;
        int lastPhase, count, phaseCount;
        int minVecX;
        int maxVecX;
        bool boundedVecX;

        DerivedFunc(const T &e) : 
            expr(e), lastPhase(-1), count(-1), phaseCount(0),
            minVecX(e.minVecX()), maxVecX(e.maxVecX()), boundedVecX(e.boundedVecX()) {
            size[0] = e.getSize(0);
            size[1] = e.getSize(1);
            size[2] = e.getSize(2);
            size[3] = e.getSize(3);
        }

        void realize(Image m) {
            m.set(expr);
        }

        void evalScanline(int y, int t, int c) {
            
            //printf("Evaluating %s(%p) at scanline %d-%d %d %d %d\n",
            //name.c_str(), this, minX, maxX, y, t, c);

            /*

            assert(im.defined(), 
                   "I was not prepared for this - no memory is allocated!");

            assert(y >= minY && y < maxY && 
                   t >= minT && t < maxT && 
                   c >= minC && c < maxC, 
                   "I was not prepared for evaluation over this domain!\n");
            */

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

        void prepare(Region r, int phase) {
            /*printf("Preparing %s(%p) phase %d over region %d %d %d %d  %d %d %d %d\n", 
                   name.c_str(), this, phase, 
                   r.x, r.y, r.t, r.c, r.width, r.height, r.frames, r.channels);               
            */

            if (phase == 0) {
                // Topology discovery: How many times will prepare be called for each phase
                if (lastPhase != 0) {
                    count = 1;
                    expr.prepare(r, phase);
                } else {
                    count++;
                }
                phaseCount = count;
            } else if (phase == 1) {
                // Region tracking and allocation
                if (lastPhase != 1) {
                    phaseCount = 1;
                    minX = r.x;
                    minY = r.y;
                    minT = r.t;
                    minC = r.c;
                    maxX = r.x + r.width;
                    maxY = r.y + r.height;
                    maxT = r.t + r.frames;
                    maxC = r.c + r.channels;
                } else {
                    phaseCount++;
                    minX = std::min(r.x, minX);
                    minY = std::min(r.y, minY);
                    minT = std::min(r.t, minT);
                    minC = std::min(r.c, minC);
                    maxX = std::max(r.x + r.width, maxX);
                    maxY = std::max(r.y + r.height, maxY);
                    maxT = std::max(r.t + r.frames, maxT);
                    maxC = std::max(r.c + r.channels, maxC);
                }
                if (phaseCount == count) {
                    if (!im.defined() ||                            
                        im.width < maxX - minX ||
                        im.height < maxY - minY ||
                        im.frames < maxT - minT ||
                        im.channels < maxC - minC) {
                        //printf("Allocating backing for %p %d %d %d %d\n", 
                        //this, maxX-minX, maxY-minY, maxT-minT, maxC-minC);
                        // We over-allocate by Vec::width so that vectorizing is always safe (boundedVecX() is false)
                        im = Image(maxX - minX + Vec::width, maxY - minY, maxT - minT, maxC - minC);
                        //printf("Done\n");
                        
                    } else {
                        // no need to allocate
                    }

                    r.x = minX;
                    r.y = minY;
                    r.t = minT;
                    r.c = minC;
                    r.width = maxX - minX;
                    r.height = maxY - minY;
                    r.frames = maxT - minT;
                    r.channels = maxC - minC;
                    expr.prepare(r, 1);

                }
            } else if (phase == 2) {
                // Evaluation
                if (lastPhase != 2) {
                    phaseCount = 1;
                    expr.prepare(r, 2);

                    if (lazy) {
                        // No scanlines have been evaluated
                        evaluated.assign((maxY - minY)*(maxT - minT)*(maxC - minC), false);
                    } else {
                        /*
                        printf("Evaluating %s(%p) over %d %d %d %d %d %d %d %d\n",
                               name.c_str(), this, minX, minY, minT, minC,
                               maxX-minX, maxY-minY, maxT-minT, maxC-minC);
                        */
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
                        //printf("Done evaluating %s(%p)\n", name.c_str(), this);
                    }                                
                } else {
                    phaseCount++;
                }

            } else if (phase == 3) {
                // clean up
                im = Image();
            }                       

            lastPhase = phase;

        }

        std::pair<float, float> bounds(Expr::Region r) const {
            return expr.bounds(r);
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
        class FuncRef {
            // We shift the args to index into the image correctly,
            // because the backing for a function may have an
            // arbitrary offset
            const std::shared_ptr<BaseFunc> ptr;
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
            }

            int getSize(int i) const {
                return sizes[i];
            }

            typedef ImRefIter<IBinaryOp<SX, ConstInt, Vec::Sub>, 
                              IBinaryOp<SY, ConstInt, Vec::Sub>, 
                              IBinaryOp<ST, ConstInt, Vec::Sub>, 
                              IBinaryOp<SC, ConstInt, Vec::Sub>, 
                              AffineCase, ShiftedCase> Iter;

            Iter scanline(int x, int y, int t, int c, int width) const {
                auto sxIter = (sx-ptr->minX).scanline(x, y, t, c, width);
                auto syIter = (sy-ptr->minY).scanline(x, y, t, c, width);
                auto stIter = (st-ptr->minT).scanline(x, y, t, c, width);
                auto scIter = (sc-ptr->minC).scanline(x, y, t, c, width);
                if (ptr->lazy) {
                    if (AffineCase || ShiftedCase) {
                        ptr->evalScanlineIfNeeded(syIter[0]+ptr->minY, 
                                                  stIter[0]+ptr->minT, 
                                                  scIter[0]+ptr->minC);                    
                    } else {
                        Region r = {x, y, t, c, width, 1, 1, 1};
                        std::pair<int, int> yb = sy.bounds(r);
                        std::pair<int, int> tb = st.bounds(r);
                        std::pair<int, int> cb = sc.bounds(r);
                        /*
                          printf("Evaluating over %d %d %d  %d %d %d\n", 
                          yb.first, tb.first, cb.first, 
                          yb.second, tb.second, cb.second);
                        */
                        for (int ec = cb.first; ec <= cb.second; ec++) {
                            for (int et = tb.first; et <= tb.second; et++) {
                                for (int ey = yb.first; ey <= yb.second; ey++) {
                                    ptr->evalScanlineIfNeeded(ey, et, ec);
                                }
                            }
                        }
                    }
                }
                return Iter(ptr->im, sxIter, syIter, stIter, scIter);
            }


            // We over-allocate image backing so that this is always safe
            bool boundedVecX() const {
                return false;
            }
            int minVecX() const {
                return -HUGE_INT;
            }
            int maxVecX() const {
                return HUGE_INT;
            }

            std::pair<float, float> bounds(Region r) const {
                //return ptr->bounds(r); 

                // we could recurse here, but we intentionally don't
                // to prevent an exponential explosion of bounds calls
                return std::make_pair(-INF, INF);
            }

            void prepare(Region r, int phase) const {
                // We require whatever the args reference
                sx.prepare(r, phase);
                sy.prepare(r, phase);
                st.prepare(r, phase);
                sc.prepare(r, phase);

                // Plus the function itself over the bounds of what the args could evaluate to
                std::pair<int, int> xb = sx.bounds(r);
                std::pair<int, int> yb = sy.bounds(r);
                std::pair<int, int> tb = st.bounds(r);
                std::pair<int, int> cb = sc.bounds(r);

                Region r2 = {xb.first, yb.first, tb.first, cb.first, 
                             xb.second - xb.first + 1,
                             yb.second - yb.first + 1,
                             tb.second - tb.first + 1,
                             cb.second - cb.first + 1};
                ptr->prepare(r2, phase);
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
            return ptr->size[i];
        }

        bool boundedVecX() const {return false;}
        int minVecX() const {return -HUGE_INT;} 
        int maxVecX() const {return HUGE_INT;}

        void prepare(Region r, int phase) const {
            ptr->prepare(r, phase);
        }

        std::pair<float, float> bounds(Region r) const {
            //return ptr->bounds(r);
            return std::make_pair(-INF, INF);
        }
        
        void lazy() {
            ptr->lazy = true;
        }
        
        void eager() {
            ptr->lazy = false;
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

        // For funcs with bounds, evaluate over the bounds into a new image
        Image realize() const {
            assert(getSize(0) && getSize(1) && getSize(2) && getSize(3), 
                   "Can only construct an image from a bounded function.\n"
                   "Call realize(width, height, frames, channels) instead.\n");
            Image im(getSize(0), getSize(1), getSize(2), getSize(3));
            ptr->realize(im);
            return im;
        }

        // For debugging
        void setName(const std::string &str) {
            ptr->name = str;
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
