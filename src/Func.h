#ifndef IMAGESTACK_FUNC_H
#define IMAGESTACK_FUNC_H

#include "Image.h"
#include "header.h"

namespace Lazy {

    struct BaseFunc {
        Image m;
        int xo, yo, to, co;

        vector<pair<int, int> > evaluated;

        // Evaluate myself into buffer at the given offset. 
        virtual void evalScanline(int x, int y, int t, int c, int width) = 0;

        // Prepare to be evaluated over a given region
        virtual void prepareFunc(int x, int y, int t, int c, 
                                        int width, int height, int frames, int channels) = 0;

        virtual int getSize(int i) = 0;
    };

    template<typename T>
    struct DerivedFunc : public BaseFunc {
        const T expr;
    
        DerivedFunc(const T e) : expr(e) {}

        void evalScanline(int x, int y, int t, int c, int width) {
            // TODO: reuse?
            printf("Evaluating at scanline %d %d %d %d %d\n", 
                   x, y, t, c, width);        
            m.setScanline(expr.scanline(x, y, t, c, width), x-xo, y-yo, t-to, c-co, width);
        }   

        void prepareFunc(int x, int y, int t, int c,
                                int width, int height, int frames, int channels) {
            printf("Preparing over region %d %d %d %d  %d %d %d %d\n",
                   x, y, t, c, width, height, frames, channels);               

            xo = x; yo = y; to = t; co = c;

            // TODO: doesn't work for multiple shifted usages. Only
            // allocates enough for the last usage. Four-phases? (clear, accumulate, finalize, free)
            if (m.defined() &&
                m.width == width &&
                m.height == height &&
                m.frames == frames &&
                m.channels == channels) {
                // We're good - already allocated
            } else {
                m = Image(width, height, frames, channels);
                evaluated.resize(height*frames*channels);
                memset(&evaluated[0], 0, evaluated.size() * sizeof(evaluated[0]));
            }
            expr.prepare(x, y, t, c, width, height, frames, channels);
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
        }

        // Implement the lazy interface
        typedef Func Lazy;
        float operator()(int x, int y, int t, int c) const {
            return ptr->m(x - ptr->xo, y - ptr->yo, t - ptr->to , c - ptr->co);
        }

        int getSize(int i) const {
            return ptr->getSize(i);
        }

        typedef Image::Iter Iter;
        Iter scanline(int x, int y, int t, int c, int width) const {        
            pair<int, int> &run = ptr->evaluated[(c*ptr->m.frames + t)*ptr->m.height + y];

            if (run.second == 0) {
                // There is nothing already evaluated
                ptr->evalScanline(x, y, t, c, width);
                run.first = x;
                run.second = x+width;
            } else {
                bool startsAfter = x >= run.first;
                bool endsBefore = x + width <= run.second;

                // Compute more stuff at the start
                if (x < run.first) {
                    ptr->evalScanline(x, y, t, c, run.first - x);
                    run.first = x;
                }

                // Compute more stuff at the end
                if (x + width > run.second) {
                    ptr->evalScanline(run.second, y, t, c, x + width - run.second);
                    run.second = x + width;
                }
            }
            return ptr->m.scanline(x-ptr->xo, y-ptr->yo, t-ptr->yo, c-ptr->co, width);
        }
    
        void prepare(int x, int y, int t, int c,
                     int width, int height, int frames, int channels) const {
            ptr->prepareFunc(x, y, t, c,
                             width, height, frames, channels);
        }
    };
}

#include "footer.h"
#endif
