#ifndef IMAGESTACK_FUNC_H
#define IMAGESTACK_FUNC_H

#include "Image.h"
#include "header.h"

namespace Lazy {

    struct BaseFunc {
        Image im;
        int minX, minY, minT, minC;
        int maxX, maxY, maxT, maxC;
        bool lazy;

        vector<int> evaluated;

        // Evaluate myself into buffer at the given scanline. 
        virtual void evalScanline(int y, int t, int c) = 0;

        // Prepare to be evaluated over a given region
        virtual void prepareFunc(int phase, int x, int y, int t, int c, 
                                 int width, int height, int frames, int channels) = 0;

        virtual int getSize(int i) = 0;
    };

    template<typename T>
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
            //printf("Evaluating %p at scanline %d-%d %d %d %d (placing in image at %d-%d %d %d %d)\n", 
            // this, minX, maxX, y, t, c, 0, maxX-minX, y-minY, t-minT, c-minC);        
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

        void prepareFunc(int phase, int x, int y, int t, int c,
                         int width, int height, int frames, int channels) {

            // recurse
            expr.prepare(phase, x, y, t, c, width, height, frames, channels);

            //printf("Preparing %p at phase %d over region %d %d %d %d  %d %d %d %d\n", this,
            //phase, x, y, t, c, width, height, frames, channels);               

            // Each phase we get called multiple times according to
            // how many times we occur in the expression

            // The start of a new page
            // In phase zero we take unions of the regions requested
            if (phase == 0) {
                if (lastPhase != 0) {
                    minX = x;
                    minY = y;
                    minT = t;
                    minC = c;
                    maxX = x + width;
                    maxY = y + height;
                    maxT = t + frames;
                    maxC = c + channels;
                } else {
                    minX = std::min(minX, x);
                    minY = std::min(minY, y);
                    minT = std::min(minT, t);
                    minC = std::min(minC, c);
                    maxX = std::max(maxX, x + width);
                    maxY = std::max(maxY, y + height);
                    maxT = std::max(maxT, t + frames);
                    maxC = std::max(maxC, c + channels);
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
                        for (int c = minC; c < maxC; c++) {
                            for (int t = minT; t < maxT; t++) {
                                #ifdef _OPENMP
                                #pragma omp parallel for
                                #endif
                                for (int y = minY; y < maxY; y++) {
                                    evalScanline(y, t, c);
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
        Func(const T t, const typename T::Lazy *enable = NULL) {
            ptr.reset(new DerivedFunc<T>(t));
            ptr->lazy = true;
        }

        // Implement the lazy interface
        typedef Func Lazy;

        int getSize(int i) const {
            return ptr->getSize(i);
        }

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

        void prepare(int phase, int x, int y, int t, int c,
                     int width, int height, int frames, int channels) const {
            ptr->prepareFunc(phase, x, y, t, c,
                             width, height, frames, channels);
        }
    };
}

#include "footer.h"
#endif
