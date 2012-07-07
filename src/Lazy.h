#ifndef IMAGESTACK_LAZY_H
#define IMAGESTACK_LAZY_H

#include <immintrin.h>
#include <stdint.h>
#include "header.h"

// This file defines a set of image-like function objects. They
// represent pure functions over a 4-d integer domain. Their body is
// completely contained within their type, so that they compile to
// efficient code (with no dynamic dispatch).
//
// E.g., the function f(x, y, t, c) = x*3 + 4 has type Add<Mul<X, Const>, Const>
//
// They are all tagged with a nested type called Lazy so that sfinae prevents unwanted constructions

class Image;
namespace Lazy {

    namespace Vec {
#ifdef __AVX__
        typedef __m256 type;
        const int width = 8;

        static type broadcast(float v) {
            return _mm256_set1_ps(v);
        }

        static type zero() {
            return _mm256_setzero_ps();
        }

        // Arithmetic binary operators
        struct Add {
            static float scalar(float a, float b) {return a + b;}
            static type vec(type a, type b) {return _mm256_add_ps(a, b);}
        };
        struct Sub {
            static float scalar(float a, float b) {return a - b;}
            static type vec(type a, type b) {return _mm256_sub_ps(a, b);}
        };
        struct Mul {
            static float scalar(float a, float b) {return a * b;}
            static type vec(type a, type b) {return _mm256_mul_ps(a, b);}
        };
        struct Div {
            static float scalar(float a, float b) {return a / b;}
            static type vec(type a, type b) {return _mm256_div_ps(a, b);}
        };
        struct Min {
            static float scalar(float a, float b) {return a < b ? a : b;}
            static type vec(type a, type b) {return _mm256_min_ps(a, b);}
        };
        struct Max {
            static float scalar(float a, float b) {return a > b ? a : b;}
            static type vec(type a, type b) {return _mm256_max_ps(a, b);}
        };

        // Comparisons
        struct GT {
            static bool scalar(float a, float b) {return a > b;}
            static type vec(type a, type b) {return _mm256_cmp_ps(a, b, _CMP_GT_OQ);}
        };
        struct LT {
            static bool scalar(float a, float b) {return a < b;}
            static type vec(type a, type b) {return _mm256_cmp_ps(a, b, _CMP_LT_OQ);}
        };
        struct GE {
            static bool scalar(float a, float b) {return a >= b;}
            static type vec(type a, type b) {return _mm256_cmp_ps(a, b, _CMP_GE_OQ);}
        };
        struct LE {
            static bool scalar(float a, float b) {return a <= b;}
            static type vec(type a, type b) {return _mm256_cmp_ps(a, b, _CMP_LE_OQ);}
        };
        struct EQ {
            static bool scalar(float a, float b) {return a == b;}
            static type vec(type a, type b) {return _mm256_cmp_ps(a, b, _CMP_EQ_OQ);}
        };
        struct NEQ {
            static bool scalar(float a, float b) {return a != b;}
            static type vec(type a, type b) {return _mm256_cmp_ps(a, b, _CMP_NEQ_OQ);}
        };

        // Logical ops
        static type blend(type a, type b, type mask) {
            return _mm256_blendv_ps(a, b, mask);
        }

        static type interleave(type a, type b) {
            __m256 r_lo = _mm256_unpacklo_ps(a, b);
            __m256 r_hi = _mm256_unpackhi_ps(a, b);
            return _mm256_permute2f128_ps(r_lo, r_hi, 1<<5);
        }

        // Unary ops
        struct Floor {
            static float scalar(float a) {return floorf(a);}
            static type vec(type a) {return _mm256_floor_ps(a);}
        };
        struct Ceil {
            static float scalar(float a) {return ceilf(a);}
            static type vec(type a) {return _mm256_ceil_ps(a);}
        };
        struct Sqrt {
            static float scalar(float a) {return sqrtf(a);}
            static type vec(type a) {return _mm256_sqrt_ps(a);}
        };

        // Loads and stores
        static type load(const float *f) {
            return _mm256_loadu_ps(f);
        }

        /*
        static type maskedLoad(const float *f, type mask) {
            return _mm256_maskload_ps(f, mask);
        }
        */

        static void store(type a, float *f) {
            _mm256_storeu_ps(f, a);
        }

#else
#ifdef __SSE__
        typedef __m128 type;
        const int width = 4;

        static type broadcast(float v) {
            return _mm_set1_ps(v);
        }

        static type zero() {
            return _mm_setzero_ps();
        }

        // Arithmetic binary operators
        struct Add {
            static float scalar(float a, float b) {return a + b;}
            static type vec(type a, type b) {return _mm_add_ps(a, b);}
        };
        struct Sub {
            static float scalar(float a, float b) {return a - b;}
            static type vec(type a, type b) {return _mm_sub_ps(a, b);}
        };
        struct Mul {
            static float scalar(float a, float b) {return a * b;}
            static type vec(type a, type b) {return _mm_mul_ps(a, b);}
        };
        struct Div {
            static float scalar(float a, float b) {return a / b;}
            static type vec(type a, type b) {return _mm_div_ps(a, b);}
        };
        struct Min {
            static float scalar(float a, float b) {return a < b ? a : b;}
            static type vec(type a, type b) {return _mm_min_ps(a, b);}
        };
        struct Max {
            static float scalar(float a, float b) {return a > b ? a : b;}
            static type vec(type a, type b) {return _mm_max_ps(a, b);}
        };

        // Comparisons
        struct GT {
            static bool scalar(float a, float b) {return a > b;}
            static type vec(type a, type b) {return _mm_cmpgt_ps(a, b);}
        };
        struct LT {
            static bool scalar(float a, float b) {return a < b;}
            static type vec(type a, type b) {return _mm_cmplt_ps(a, b);}
        };
        struct GE {
            static bool scalar(float a, float b) {return a >= b;}
            static type vec(type a, type b) {return _mm_cmpge_ps(a, b);}
        };
        struct LE {
            static bool scalar(float a, float b) {return a <= b;}
            static type vec(type a, type b) {return _mm_cmple_ps(a, b);}
        };
        struct EQ {
            static bool scalar(float a, float b) {return a == b;}
            static type vec(type a, type b) {return _mm_cmpeq_ps(a, b);}
        };
        struct NEQ {
            static bool scalar(float a, float b) {return a != b;}
            static type vec(type a, type b) {return _mm_cmpneq_ps(a, b);}
        };

        static type interleave(type a, type b) {
            return _mm_unpacklo_ps(a, b);
        }

#ifdef __SSE4_1__
        // Logical ops
        static type blend(type a, type b, type mask) {
            return _mm_blendv_ps(a, b, mask);
        }

        // Unary ops
        struct Floor {
            static float scalar(float a) {return floorf(a);}
            static type vec(type a) {return _mm_floor_ps(a);}
        };
        struct Ceil {
            static float scalar(float a) {return ceilf(a);}
            static type vec(type a) {return _mm_ceil_ps(a);}
        };
        struct Sqrt {
            static float scalar(float a) {return sqrtf(a);}
            static type vec(type a) {return _mm_sqrt_ps(a);}
        };
#else

        static type blend(type a, type b, type mask) {
            return _mm_or_ps(_mm_and_ps(mask, b),
                             _mm_andnot_ps(mask, a));
        }

        struct Floor {
            static float scalar(float a) {return floorf(a);}
            static type vec(type a) {
                union {
                    float f[width];
                    type v;
                } v;
                v.v = a;
                v.f[0] = scalar(v.f[0]);
                v.f[1] = scalar(v.f[1]);
                v.f[2] = scalar(v.f[2]);
                v.f[3] = scalar(v.f[3]);
                return v.v;
            }
        };
        struct Ceil {
            static float scalar(float a) {return ceilf(a);}
            static type vec(type a) {
                union {
                    float f[width];
                    type v;
                } v;
                v.v = a;
                v.f[0] = scalar(v.f[0]);
                v.f[1] = scalar(v.f[1]);
                v.f[2] = scalar(v.f[2]);
                v.f[3] = scalar(v.f[3]);
                return v.v;
            }
        };
        struct Sqrt {
            static float scalar(float a) {return sqrtf(a);}
            static type vec(type a) {
                union {
                    float f[width];
                    type v;
                } v;
                v.v = a;
                v.f[0] = scalar(v.f[0]);
                v.f[1] = scalar(v.f[1]);
                v.f[2] = scalar(v.f[2]);
                v.f[3] = scalar(v.f[3]);
                return v.v;
            }
        };

#endif

        // Loads and stores
        static type load(const float *f) {
            return _mm_loadu_ps(f);
        }

        static void store(type a, float *f) {
            _mm_storeu_ps(f, a);
        }

#else
        // scalar fallback
        typedef float type;
        const int width = 1;

        static type broadcast(float v) {
            return v;
        }

        static type zero() {
            return 0;
        }

        // Arithmetic binary operators
        struct Add {
            static float scalar(float a, float b) {return a + b;}
            static type vec(type a, type b) {return scalar(a, b);}
        };
        struct Sub {
            static float scalar(float a, float b) {return a - b;}
            static type vec(type a, type b) {return scalar(a, b);}
        };
        struct Mul {
            static float scalar(float a, float b) {return a * b;}
            static type vec(type a, type b) {return scalar(a, b);}
        };
        struct Div {
            static float scalar(float a, float b) {return a / b;}
            static type vec(type a, type b) {return scalar(a, b);}
        };
        struct Min {
            static float scalar(float a, float b) {return a < b ? a : b;}
            static type vec(type a, type b) {return scalar(a, b);}
        };
        struct Max {
            static float scalar(float a, float b) {return a > b ? a : b;}
            static type vec(type a, type b) {return scalar(a, b);}
        };

        // Comparisons
        struct GT {
            static bool scalar(float a, float b) {return a > b;}
            static bool vec(type a, type b) {return scalar(a, b);}
        };
        struct LT {
            static bool scalar(float a, float b) {return a < b;}
            static bool vec(type a, type b) {return scalar(a, b);}
        };
        struct GE {
            static bool scalar(float a, float b) {return a >= b;}
            static bool vec(type a, type b) {return scalar(a, b);}
        };
        struct LE {
            static bool scalar(float a, float b) {return a <= b;}
            static bool vec(type a, type b) {return scalar(a, b);}
        };
        struct EQ {
            static bool scalar(float a, float b) {return a == b;}
            static bool vec(type a, type b) {return scalar(a, b);}
        };
        struct NEQ {
            static bool scalar(float a, float b) {return a != b;}
            static bool vec(type a, type b) {return scalar(a, b);}
        };

        // Logical ops
        static type blend(bool mask, type a, type b) {
            return (mask ? b : a);
        }

        static type interleave(type a, type b) {
            return a;
        }

        // Unary ops
        struct Floor {
            static float scalar(float a) {return floorf(a);}
            static type vec(type a) {return scalar(a);}
        };
        struct Ceil {
            static float scalar(float a) {return ceilf(a);}
            static type vec(type a) {return scalar(a);}
        };
        struct Sqrt {
            static float scalar(float a) {return sqrtf(a);}
            static type vec(type a) {return scalar(a);}
        };

        // Loads and stores
        static type load(const float *f) {
            return *f;
        }

        static void store(type a, float *f) {
            *f = a;
        }


#endif
#endif
    }

    // How should expressions hold onto sub-expressions? We use
    // const copies. Experiments were made with const references to
    // try and trick g++ into doing more aggressive subexpression
    // elimination, but results were mixed at best.
    template<typename T>
    struct Handle {
        typedef const T type;
    };

    // A base class for things which do not depend on image data
    struct Unbounded {
        int getSize(int) const {return 0;}
    };

    // Constants
    struct Const : public Unbounded {
        typedef Const Lazy;
        const float val;
        Const(const float val_) : val(val_) {}
        float operator()(int x, int, int, int) const {return val;}

        // State needed to iterate across a scanline
        struct Iter {
            const float val;
            const Vec::type vec_val;
            Iter() : val(0), vec_val(Vec::zero()) {}
            Iter(float v) : val(v), vec_val(Vec::broadcast(val)) {
            }
            float operator[](int x) const {return val;}
            Vec::type vec(int x) const {return vec_val;}
        };
        Iter scanline(int x, int y, int t, int c, int width) const {
            return Iter(val);
        }
        void prepare(int x, int y, int t, int c, 
                     int width, int height, int frames, int channels) const {}
    };

    template<>
    struct Handle<Const> {
        typedef const Const type;
    };

    // Coordinates
    struct X : public Unbounded {
        typedef X Lazy;
        float operator()(int x, int, int, int) const {return x;}

        // State needed to iterate across a scanline
        struct Iter {
            float operator[](int x) const {return x;}
            Vec::type vec(int x) const {
                union {
                    float f[Vec::width];
                    Vec::type v;
                } v;
                for (int i = 0; i < Vec::width; i++) {
                    v.f[i] = x+i;
                }
                return v.v;
            }
        };
        Iter scanline(int x, int y, int t, int c, int width) const {
            return Iter();
        }
        void prepare(int x, int y, int t, int c, 
                     int width, int height, int frames, int channels) const {}
    };

    struct Y : public Unbounded {
        typedef Y Lazy;
        float operator()(int, int y, int, int) const {return y;}

        typedef Const::Iter Iter;

        Const::Iter scanline(int x, int y, int t, int c, int width) const {
            return Const::Iter(y);
        }
        void prepare(int x, int y, int t, int c, 
                     int width, int height, int frames, int channels) const {}
    };

    struct T : public Unbounded {
        typedef T Lazy;
        float operator()(int, int, int t, int) const {return t;}

        typedef Const::Iter Iter;

        Const::Iter scanline(int x, int y, int t, int c, int width) const {
            return Const::Iter(t);
        }
        void prepare(int x, int y, int t, int c, 
                     int width, int height, int frames, int channels) const {}
    };

    struct C : public Unbounded {
        typedef C Lazy;
        float operator()(int, int, int, int c) const {return c;}

        typedef Const::Iter Iter;

        Const::Iter scanline(int x, int y, int t, int c, int width) const {
            return Const::Iter(c);
        }
        void prepare(int x, int y, int t, int c, 
                     int width, int height, int frames, int channels) const {}
    };

    // Arithmetic binary operators
    template<typename A, typename B, typename Op>
    struct BinaryOp {
        typedef BinaryOp<typename A::Lazy, typename B::Lazy, Op> Lazy;
        typename Handle<A>::type a;
        typename Handle<B>::type b;
        
        BinaryOp(const A &a_, const B &b_) : a(a_), b(b_) {
            for (int i = 0; i < 4; i++) {
                if (a.getSize(i) && b.getSize(i)) {
                    assert(a.getSize(i) == b.getSize(i),
                           "Can only combine images with matching size\n");
                }
            }
        }
        
        int getSize(int i) const {
            if (a.getSize(i)) return a.getSize(i);
            return b.getSize(i);
        }
        
        float operator()(int x, int y, int t, int c) const {
            return Op::scalar(a(x, y, t, c), b(x, y, t, c));
        }
        
        struct Iter {
            const typename A::Iter a;
            const typename B::Iter b;
            Iter() {}
            Iter(const typename A::Iter &a_, const typename B::Iter &b_) : a(a_), b(b_) {}
            float operator[](int x) const {
                return Op::scalar(a[x], b[x]);
            }
            Vec::type vec(int x) const {
                return Op::vec(a.vec(x), b.vec(x));
            }
        };
        Iter scanline(int x, int y, int t, int c, int width) const {
            return Iter(a.scanline(x, y, t, c, width), b.scanline(x, y, t, c, width));
        }
        void prepare(int x, int y, int t, int c, 
                     int width, int height, int frames, int channels) const {
            a.prepare(x, y, t, c, 
                      width, height, frames, channels);
            b.prepare(x, y, t, c, 
                      width, height, frames, channels);
        }
    };
    
    // Comparison binary operators
    template<typename A, typename B, typename Op>
    struct Cmp {
        typedef Cmp<typename A::Lazy, typename B::Lazy, Op> LazyBool;
        typename Handle<A>::type a;
        typename Handle<B>::type b;
        
        Cmp(const A &a_, const B &b_) : a(a_), b(b_) {
            for (int i = 0; i < 4; i++) {
                if (a.getSize(i) && b.getSize(i)) {
                    assert(a.getSize(i) == b.getSize(i),
                           "Can only combine images with matching size\n");
                }
            }
        }
        int getSize(int i) const {
            if (a.getSize(i)) return a.getSize(i);
            return b.getSize(i);
        }
        
        float operator()(int x, int y, int t, int c) const {
            return Op::scalar(a(x, y, t, c), b(x, y, t, c)) ? 1 : 0;
        }
        
        struct Iter {
            const typename A::Iter a;
            const typename B::Iter b;
            Iter() {}
            Iter(const typename A::Iter &a_, const typename B::Iter &b_) : a(a_), b(b_) {}
            float operator[](int x) const {
                return Op::scalar(a[x], b[x]) ? 1 : 0;
            }
            Vec::type vec(int x) const {
                return Op::vec(a.vec(x), b.vec(x));
            }
        };
        Iter scanline(int x, int y, int t, int c, int width) const {
            return Iter(a.scanline(x, y, t, c, width),
                        b.scanline(x, y, t, c, width));
        }
        void prepare(int x, int y, int t, int c, 
                     int width, int height, int frames, int channels) const {
            a.prepare(x, y, t, c, 
                      width, height, frames, channels);
            b.prepare(x, y, t, c, 
                      width, height, frames, channels);
        }
    };
    
    template<typename A, typename B>
    BinaryOp<typename A::Lazy, typename B::Lazy, Vec::Min>
    min(const A &a, const B &b) {
        return BinaryOp<typename A::Lazy, typename B::Lazy, Vec::Min>(a, b);
    }
    template<typename A>
    BinaryOp<typename A::Lazy, Const, Vec::Min>
    min(const A &a, float b) {
        return min(a, Const(b));
    }
    template<typename B>
    BinaryOp<Const, typename B::Lazy, Vec::Min>
    min(float a, const B &b) {
        return min(Const(a), b);
    }
    template<typename A, typename B>
    BinaryOp<typename A::Lazy, typename B::Lazy, Vec::Max>
    max(const A &a, const B &b) {
        return BinaryOp<typename A::Lazy, typename B::Lazy, Vec::Max>(a, b);
    }
    template<typename A>
    BinaryOp<typename A::Lazy, Const, Vec::Max>
    max(const A &a, float b) {
        return max(a, Const(b));
    }
    template<typename B>
    BinaryOp<Const, typename B::Lazy, Vec::Max>
    max(float a, const B &b) {
        return max(Const(a), b);
    }

    template<typename A, typename B, typename C>
    BinaryOp<BinaryOp<typename A::Lazy, typename B::Lazy, Vec::Max>, typename C::Lazy, Vec::Min>
    clamp(const A &a, const B &b, const C &c) {
        return min(max(a, b), c);
    }
    template<typename B, typename C>
    BinaryOp<BinaryOp<Const, typename B::Lazy, Vec::Max>, typename C::Lazy, Vec::Min>
    clamp(float a, const B &b, const C &c) {
        return min(max(a, b), c);
    }
    template<typename A, typename C>
    BinaryOp<BinaryOp<typename A::Lazy, Const, Vec::Max>, typename C::Lazy, Vec::Min>
    clamp(const A &a, float b, const C &c) {
        return min(max(a, b), c);
    }
    template<typename A, typename B>
    BinaryOp<BinaryOp<typename A::Lazy, typename B::Lazy, Vec::Max>, Const, Vec::Min>
    clamp(const A &a, const B &b, float c) {
        return min(max(a, b), c);
    }
    template<typename A>
    BinaryOp<BinaryOp<typename A::Lazy, Const, Vec::Max>, Const, Vec::Min>
    clamp(const A &a, float b, float c) {
        return min(max(a, b), c);
    }

    // Lift a unary function over floats to the same function over an image (e.g. cosf)
    template<float(*fn)(float), typename A>
    struct Lift {
        typedef Lift<fn, typename A::Lazy> Lazy;
        typename Handle<A>::type a;
        Lift(const A &a_) : a(a_) {}
        float operator()(int x, int y, int t, int c) const {
            return (*fn)(a(x, y, t, c));
        }

        int getSize(int i) const {return a.getSize(i);}

        struct Iter {
            const typename A::Iter a;
            Iter() {}
            Iter(const typename A::Iter &a_) : a(a_) {}
            float operator[](int x) const {return (*fn)(a[x]);}
            Vec::type vec(int x) const {
                union {
                    float f[Vec::width];
                    Vec::type v;
                } va;
                va.v = a.vec(x);
                for (int i = 0; i < Vec::width; i++) {
                    va.f[i] = (*fn)(va.f[i]);
                }
                return va.v;
            }
        };
        Iter scanline(int x, int y, int t, int c, int width) const {
            return Iter(a.scanline(x, y, t, c, width));
        }
        void prepare(int x, int y, int t, int c, 
                     int width, int height, int frames, int channels) const {
            a.prepare(x, y, t, c, 
                      width, height, frames, channels);
        }
    };

    // Lift a vector function to the same function over an image (e.g. floor)
    template<typename A, typename Op>
    struct UnaryOp {
        typedef UnaryOp<typename A::Lazy, Op> Lazy;
        typename Handle<A>::type a;
        UnaryOp(const A &a_) : a(a_) {}
        float operator()(int x, int y, int t, int c) const {
            return Op::scalar(x, y, t, c);
        }

        int getSize(int i) const {return a.getSize(i);}

        struct Iter {
            const typename A::Iter a;
            Iter() {}
            Iter(const typename A::Iter &a_) : a(a_) {}
            float operator[](int x) const {return Op::scalar(a[x]);}
            Vec::type vec(int x) const {
                return Op::vec(a.vec(x));
            }
        };
        Iter scanline(int x, int y, int t, int c, int width) const {
            return Iter(a.scanline(x, y, t, c, width));
        }
        void prepare(int x, int y, int t, int c, 
                     int width, int height, int frames, int channels) const {
            a.prepare(x, y, t, c, 
                      width, height, frames, channels);
        }
    };

    // Arithmetic binary operators
    template<float(*fn)(float, float), typename A, typename B>
    struct Lift2 {
        typedef Lift2<fn, typename A::Lazy, typename B::Lazy> Lazy;
        typename Handle<A>::type a;
        typename Handle<B>::type b;

        Lift2(const A &a_, const B &b_) : a(a_), b(b_) {
            for (int i = 0; i < 4; i++) {
                if (a.getSize(i) && b.getSize(i)) {
                    assert(a.getSize(i) == b.getSize(i),
                           "Can only combine images with matching size\n");
                }
            }
        }

        int getSize(int i) const {
            if (a.getSize(i)) return a.getSize(i);
            return b.getSize(i);
        }

        float operator()(int x, int y, int t, int c) const {
            return (*fn)(a(x, y, t, c), b(x, y, t, c));
        }

        struct Iter {
            const typename A::Iter a;
            const typename B::Iter b;
            Iter() {}
            Iter(const typename A::Iter &a_, const typename B::Iter &b_) : a(a_), b(b_) {}
            float operator[](int x) const {
                return (*fn)(a[x], b[x]);
            }
            Vec::type vec(int x) const {
                union {
                    float f[Vec::width];
                    Vec::type v;
                } va, vb;
                va.v = a.vec(x);
                vb.v = b.vec(x);
                for (int i = 0; i < Vec::width; i++) {
                    vb.f[i] = (*fn)(va.f[i], vb.f[i]);
                }
                return vb.v;
            }
        };
        Iter scanline(int x, int y, int t, int c, int width) const {
            return Iter(a.scanline(x, y, t, c, width), b.scanline(x, y, t, c, width));
        }
        void prepare(int x, int y, int t, int c, 
                     int width, int height, int frames, int channels) const {
            a.prepare(x, y, t, c, 
                      width, height, frames, channels);
            b.prepare(x, y, t, c, 
                      width, height, frames, channels);
        }
    };

    template<typename A>
    Lift<logf, typename A::Lazy> log(const A &a) {
        return Lift<logf, typename A::Lazy>(a);
    }

    template<typename A>
    Lift<expf, typename A::Lazy> exp(const A &a) {
        return Lift<expf, typename A::Lazy>(a);
    }

    template<typename A>
    Lift<cosf, typename A::Lazy> cos(const A &a) {
        return Lift<cosf, typename A::Lazy>(a);
    }

    template<typename A>
    Lift<sinf, typename A::Lazy> sin(const A &a) {
        return Lift<sinf, typename A::Lazy>(a);
    }

    template<typename A>
    Lift<tanf, typename A::Lazy> tan(const A &a) {
        return Lift<tanf, typename A::Lazy>(a);
    }

    template<typename A>
    Lift<acosf, typename A::Lazy> acos(const A &a) {
        return Lift<acosf, typename A::Lazy>(a);
    }

    template<typename A>
    Lift<asinf, typename A::Lazy> asin(const A &a) {
        return Lift<asinf, typename A::Lazy>(a);
    }

    template<typename A>
    Lift<atanf, typename A::Lazy> atan(const A &a) {
        return Lift<atanf, typename A::Lazy>(a);
    }

    template<typename A>
    Lift<fabsf, typename A::Lazy> abs(const A &a) {
        return Lift<fabsf, typename A::Lazy>(a);
    }

    template<typename A>
    UnaryOp<typename A::Lazy, Vec::Sqrt> sqrt(const A &a) {
        return UnaryOp<typename A::Lazy, Vec::Sqrt>(a);
    }

    template<typename A>
    UnaryOp<typename A::Lazy, Vec::Floor> floor(const A &a) {
        return UnaryOp<typename A::Lazy, Vec::Floor>(a);
    }

    template<typename A>
    UnaryOp<typename A::Lazy, Vec::Ceil> ceil(const A &a) {
        return UnaryOp<typename A::Lazy, Vec::Ceil>(a);
    }

    template<typename A, typename B>
    Lift2<powf, typename A::Lazy, typename B::Lazy> pow(const A &a, const B &b) {
        return Lift2<powf, typename A::Lazy, typename B::Lazy>(a, b);
    }
    template<typename A>
    Lift2<powf, typename A::Lazy, Const> pow(const A &a, float b) {
        return Lift2<powf, typename A::Lazy, Const>(a, b);
    }
    template<typename B>
    Lift2<powf, Const, typename B::Lazy> pow(float a, const B &b) {
        return Lift2<powf, Const, typename B::Lazy>(a, b);
    }

    template<typename A, typename B>
    Lift2<fmodf, typename A::Lazy, typename B::Lazy> fmod(const A &a, const B &b) {
        return Lift2<fmodf, typename A::Lazy, typename B::Lazy>(a, b);
    }
    template<typename A>
    Lift2<fmodf, typename A::Lazy, Const> fmod(const A &a, float b) {
        return Lift2<fmodf, typename A::Lazy, Const>(a, b);
    }
    template<typename B>
    Lift2<fmodf, Const, typename B::Lazy> fmod(float a, const B &b) {
        return Lift2<fmodf, Const, typename B::Lazy>(a, b);
    }

    template<typename A, typename B>
    Lift2<atan2f, typename A::Lazy, typename B::Lazy> atan2(const A &a, const B &b) {
        return Lift2<atan2f, typename A::Lazy, typename B::Lazy>(a, b);
    }

    template<typename A>
    Lift2<atan2f, typename A::Lazy, Const> atan2(const A &a, float b) {
        return Lift2<atan2f, typename A::Lazy, Const>(a, b);
    }

    template<typename B>
    Lift2<atan2f, Const, typename B::Lazy> atan2(float a, const B &b) {
        return Lift2<atan2f, Const, typename B::Lazy>(a, b);
    }


    template<typename A, typename B, typename C>
    struct _Select {
        typedef _Select<typename A::LazyBool, typename B::Lazy, typename C::Lazy> Lazy;
        typename Handle<A>::type a;
        typename Handle<B>::type b;
        typename Handle<C>::type c;

        _Select(const A &a_, const B &b_, const C &c_) : a(a_), b(b_), c(c_) {
            for (int i = 0; i < 4; i++) {
                int s = a.getSize(i);
                if (!s) s = b.getSize(i);
                if (!s) s = c.getSize(i);
                assert((a.getSize(i) == s || a.getSize(i) == 0) &&
                       (b.getSize(i) == s || b.getSize(i) == 0) &&
                       (c.getSize(i) == s || c.getSize(i) == 0),
                       "Can only combine images with matching size\n");
            }
        }

        int getSize(int i) const {
            if (a.getSize(i)) return a.getSize(i);
            if (b.getSize(i)) return b.getSize(i);
            if (c.getSize(i)) return c.getSize(i);
            return 0;
        }

        float operator()(int x, int y, int t, int c_) const {
            return (a(x, y, t, c_) ? b(x, y, t, c_) : c(x, y, t, c_));
        }

        struct Iter {
            const typename A::Iter a;
            const typename B::Iter b;
            const typename C::Iter c;
            Iter() {}
            Iter(const typename A::Iter &a_,
                 const typename B::Iter &b_,
                 const typename C::Iter &c_) : a(a_), b(b_), c(c_) {}
            float operator[](int x) const {
                return a[x] ? b[x] : c[x];
            }
            Vec::type vec(int x) const {
                const Vec::type va = a.vec(x);
                const Vec::type vb = b.vec(x);
                const Vec::type vc = c.vec(x);
                return Vec::blend(vc, vb, va);
            }
        };
        Iter scanline(int x, int y, int t, int c_, int width) const {
            return Iter(a.scanline(x, y, t, c_, width),
                        b.scanline(x, y, t, c_, width),
                        c.scanline(x, y, t, c_, width));
        }
        void prepare(int x, int y, int t, int c_, 
                     int width, int height, int frames, int channels) const {
            a.prepare(x, y, t, c_, 
                      width, height, frames, channels);
            b.prepare(x, y, t, c_, 
                      width, height, frames, channels);
            c.prepare(x, y, t, c_, 
                      width, height, frames, channels);
        }
    };

    template<typename A, typename B, typename C>
    _Select<typename A::LazyBool, typename B::Lazy, typename C::Lazy>
    Select(const A &a, const B &b, const C &c) {
        return _Select<typename A::LazyBool, typename B::Lazy, typename C::Lazy>(a, b, c);
    }

    template<typename A, typename C>
    _Select<typename A::LazyBool, Const, typename C::Lazy>
    Select(const A &a, float b, const C &c) {
        return _Select<typename A::LazyBool, Const, typename C::Lazy>(a, Const(b), c);
    }

    template<typename A, typename B>
    _Select<typename A::LazyBool, typename B::Lazy, Const>
    Select(const A &a, const B &b, float c) {
        return _Select<typename A::LazyBool, typename B::Lazy, Const>(a, b, Const(c));
    }

    template<typename A>
    _Select<typename A::LazyBool, Const, Const>
    Select(const A &a, float b, float c) {
        return _Select<typename A::LazyBool, Const, Const>(a, Const(b), Const(c));
    }

    template<typename A, typename B, typename C>
    struct _IfThenElse {
        typedef _IfThenElse<typename A::LazyBool, typename B::Lazy, typename C::Lazy> Lazy;
        typename Handle<A>::type a;
        typename Handle<B>::type b;
        typename Handle<C>::type c;

        _IfThenElse(const A &a_, const B &b_, const C &c_) : a(a_), b(b_), c(c_) {
            for (int i = 0; i < 4; i++) {
                int s = a.getSize(i);
                if (!s) s = b.getSize(i);
                if (!s) s = c.getSize(i);
                assert((a.getSize(i) == s || a.getSize(i) == 0) &&
                       (b.getSize(i) == s || b.getSize(i) == 0) &&
                       (c.getSize(i) == s || c.getSize(i) == 0),
                       "Can only combine images with matching size\n");
            }
        }

        int getSize(int i) const {
            if (a.getSize(i)) return a.getSize(i);
            if (b.getSize(i)) return b.getSize(i);
            if (c.getSize(i)) return c.getSize(i);
            return 0;
        }

        float operator()(int x, int y, int t, int c_) const {
            return (a(x, y, t, c_) ? b(x, y, t, c_) : c(x, y, t, c_));
        }

        struct Iter {
            const typename A::Iter a;
            const typename B::Iter b;
            const typename C::Iter c;
            Iter() {}
            Iter(const typename A::Iter &a_,
                 const typename B::Iter &b_,
                 const typename C::Iter &c_) : a(a_), b(b_), c(c_) {}
            float operator[](int x) const {
                return a[x] ? b[x] : c[x];
            }
            Vec::type vec(int x) const {
                union {
                    float f[Vec::width];
                    Vec::type v;
                } vres;
                for (int i = 0; i < Vec::width; i++) {
                    if (a[x+i]) vres.f[i] = b[x+i];
                    else vres.f[i] = c[x+i];
                }
                return vres.v;
            }
        };
        Iter scanline(int x, int y, int t, int c_, int width) const {
            return Iter(a.scanline(x, y, t, c_, width),
                        b.scanline(x, y, t, c_, width),
                        c.scanline(x, y, t, c_, width));
        }
        void prepare(int x, int y, int t, int c_, 
                     int width, int height, int frames, int channels) const {
            a.prepare(x, y, t, c_, 
                      width, height, frames, channels);
            b.prepare(x, y, t, c_, 
                      width, height, frames, channels);
            c.prepare(x, y, t, c_, 
                      width, height, frames, channels);
        }
    };

    template<typename A, typename B, typename C>
    _IfThenElse<typename A::LazyBool, typename B::Lazy, typename C::Lazy>
    IfThenElse(const A &a, const B &b, const C &c) {
        return _IfThenElse<typename A::LazyBool, typename B::Lazy, typename C::Lazy>(a, b, c);
    }

    template<typename A, typename C>
    _IfThenElse<typename A::LazyBool, Const, typename C::Lazy>
    IfThenElse(const A &a, const float b, const C &c) {
        return _IfThenElse<typename A::LazyBool, Const, typename C::Lazy>(a, Const(b), c);
    }

    template<typename A, typename B>
    _IfThenElse<typename A::LazyBool, typename B::Lazy, Const>
    IfThenElse(const A &a, const B &b, const float c) {
        return _IfThenElse<typename A::LazyBool, typename B::Lazy, Const>(a, b, Const(c));
    }

    template<typename A>
    _IfThenElse<typename A::LazyBool, Const, Const>
    IfThenElse(const A &a, const float b, const float c) {
        return _IfThenElse<typename A::LazyBool, Const, Const>(a, Const(b), Const(c));
    }

    // Translation (TODO: shift by static amount using int template params?)
    template<typename A>
    struct _Shift {
        typedef _Shift<typename A::Lazy> Lazy;
        const int xo, yo, to, co;
        typename Handle<A>::type a;
        _Shift(const A &a_, int xo_, int yo_, int to_ = 0, int co_ = 0) : 
            a(a_), xo(xo_), yo(yo_), to(to_), co(co_) {
            assert((xo == 0 || a.getSize(0) == 0) &&
                   (yo == 0 || a.getSize(1) == 0) &&
                   (to == 0 || a.getSize(2) == 0) &&
                   (co == 0 || a.getSize(3) == 0),
                   "Can't shift expressions in bounded dimensions");
        }
        float operator()(int x, int y, int t, int c) const {
            return a(x-xo, y-yo, t-to, c-co);
        }
        
        int getSize(int i) const {
            return a.getSize(i);
        }
        
        struct Iter {
            const typename A::Iter a;
            const int xo;
            Iter() : xo(0) {}
            Iter(const typename A::Iter &a_, int xo_) : a(a_), xo(xo_) {}
            float operator[](int x) const {
                return a[x-xo];
            }
            Vec::type vec(int x) const {
                return a.vec(x-xo);
            }
        };
        Iter scanline(int x, int y, int t, int c, int width) const {
            return Iter(a.scanline(x-xo, y-yo, t-to, c-co, width), xo);
        }
        void prepare(int x, int y, int t, int c, 
                     int width, int height, int frames, int channels) const {
            a.prepare(x-xo, y-yo, t-to, c-co, 
                      width, height, frames, channels);
        }
    };
    
    template<typename A>
    _Shift<typename A::Lazy>
    shiftX(const A &a, int xo) {
        return _Shift<A>(a, xo, 0, 0, 0);
    }

    template<typename A>
    _Shift<typename A::Lazy>
    shiftY(const A &a, int yo) {
        return _Shift<A>(a, 0, yo, 0, 0);
    }

    template<typename A>
    _Shift<typename A::Lazy>
    shift(const A &a, int xo, int yo, int to = 0, int co = 0) {
        return _Shift<A>(a, xo, yo, to, co);
    }

    // Add a boundary condition
    template<typename A>
    struct _ZeroBoundary {
        typedef _ZeroBoundary<typename A::Lazy> Lazy;
        typename Handle<A>::type a;
        _ZeroBoundary(const A &a_) : a(a_) {}
        float operator()(int x, int y, int t, int c) const {
            if (x < 0 || x >= a.getSize(0) ||
                y < 0 || y >= a.getSize(1) ||
                t < 0 || t >= a.getSize(2) ||
                c < 0 || c >= a.getSize(3)) {
                return 0;
            } else {
                return a(x, y, t, c);
            }
        }

        // Once you add a boundary condition, things are unbounded
        int getSize(int) const {return 0;}

        struct Iter {
            const typename A::Iter a;
            const bool outOfBounds;
            const int width;
            Iter() : outOfBounds(false), width(0) {}
            Iter(const typename A::Iter &a_, int w) : 
                a(a_), outOfBounds(false), width(w) {  
            }

            float operator[](int x) const {
                if (outOfBounds || x < 0 || x >= width) return 0;
                return a[x];
            }
            Vec::type vec(int x) const {
                // Completely out-of-bounds
                if (outOfBounds || (x < 1-Vec::width) || (x >= width)) {
                    return Vec::broadcast(0);
                } else if (x >= 0 && x + Vec::width <= width) {
                    // Completely in-bounds
                    return a.vec(x);
                } else {
                    // Overlapping the start and/or end                    
                    union {
                        float f[Vec::width];
                        Vec::type v;
                    } v;
                    for (int i = 0; i < Vec::width; i++) {
                        v.f[i] = (x+i < 0 || x+i >= width) ? 0 : a[x];
                    }
                    return v.v;
                }
                
            }
        };
        Iter scanline(int x, int y, int t, int c, int width) const {
            bool oob = false;
            if (a.getSize(1)) oob = oob || y < 0 || y >= a.getSize(1);
            if (a.getSize(2)) oob = oob || t < 0 || t >= a.getSize(2);
            if (a.getSize(3)) oob = oob || c < 0 || c >= a.getSize(3);
            if (oob) return Iter();
            else {
                int xEnd = x+width;
                if (a.getSize(0)) xEnd = std::min(xEnd, a.getSize(0));
                int xStart = std::max(x, 0);
                return Iter(a.scanline(xStart, y, t, c, xEnd-xStart), a.getSize(0));
            }
        }

        void prepare(int x, int y, int t, int c, 
                     int width, int height, int frames, int channels) const {
            int xEnd = x+width;
            int yEnd = y+height;
            int tEnd = t+frames;
            int cEnd = c+channels;
            if (a.getSize(0)) xEnd = std::min(xEnd, a.getSize(0));
            if (a.getSize(1)) yEnd = std::min(yEnd, a.getSize(1));
            if (a.getSize(2)) tEnd = std::min(tEnd, a.getSize(2));
            if (a.getSize(3)) cEnd = std::min(cEnd, a.getSize(3));
            int xStart = std::max(x, 0);
            int yStart = std::max(y, 0);
            int tStart = std::max(t, 0);
            int cStart = std::max(c, 0);            
            a.prepare(xStart, yStart, tStart, cStart,
                      xEnd-xStart, yEnd-yStart, tEnd-tStart, cEnd-cStart);
        }
        
    };

    template<typename A>
    _ZeroBoundary<typename A::Lazy>
    zeroBoundary(const A &a) {
        return _ZeroBoundary<A>(a);
    }

    
    // interleavings
    template<typename A, typename B>
    struct _InterleaveX {
        typedef _InterleaveX<typename A::Lazy, typename B::Lazy> Lazy;
        typename Handle<A>::type a;
        typename Handle<B>::type b;
        
        _InterleaveX(const A &a_, const B &b_) : a(a_), b(b_) {
            for (int i = 0; i < 4; i++) {
                if (a.getSize(i) && b.getSize(i)) {
                    assert(a.getSize(i) == b.getSize(i),
                           "Can only combine images with matching size\n");
                }
            }
        }

        int getSize(int i) const {
            if (i == 0) {
                return std::max(a.getSize(0), b.getSize(0))*2;
            } else {
                return std::max(a.getSize(i), b.getSize(i));
            }
        }

        float operator()(int x, int y, int t, int c) const {
            if (x & 1) return b(x/2, y, t, c);
            else return a(x/2, y, t, c);
        }

        struct Iter {
            const typename A::Iter a;
            const typename B::Iter b;
            Iter() {}
            Iter(const typename A::Iter &a_, const typename B::Iter &b_) : a(a_), b(b_) {}
            float operator[](int x) const {
                if (x & 1) return b[x/2];
                else return a[x/2];
            }
            Vec::type vec(int x) const {                
                const int y = x/2;
                Vec::type v1, v2;
                if (x & 1) {
                    v1 = b.vec(y);
                    v2 = a.vec(y+1);
                } else {
                    v1 = a.vec(y);
                    v2 = b.vec(y);
                }
                return Vec::interleave(v1, v2);
            }
        };
 
        Iter scanline(int x, int y, int t, int c, int width) const {
            return Iter(a.scanline((x+1)/2, y, t, c, width/2), 
                        b.scanline(x/2, y, t, c, width/2));
        }

        void prepare(int x, int y, int t, int c, 
                     int width, int height, int frames, int channels) const {
            a.prepare((x+1)/2, y, t, c, width/2, height, frames, channels);
            b.prepare(x/2, y, t, c, width/2, height, frames, channels);
        }
    };

        
    template<typename A, typename B>
    _InterleaveX<typename A::Lazy, typename B::Lazy>
    interleaveX(const A &a, const B &b) {
        return _InterleaveX<A, B>(a, b);
    }

    template<typename B>
    _InterleaveX<Const, typename B::Lazy>
    interleaveX(float a, const B &b) {
        return _InterleaveX<Const, B>(Const(a), b);
    }

    template<typename A>
    _InterleaveX<typename A::Lazy, Const>
    interleaveX(const A &a, float b) {
        return _InterleaveX<A, Const>(a, Const(b));
    }

    inline _InterleaveX<Const, Const>
    interleaveX(float a, float b) {
        return _InterleaveX<Const, Const>(Const(a), Const(b));
    }

    template<typename A, typename B>
    struct _InterleaveY {
        typedef _InterleaveY<typename A::Lazy, typename B::Lazy> Lazy;
        typename Handle<A>::type a;
        typename Handle<B>::type b;
        
        _InterleaveY(const A &a_, const B &b_) : a(a_), b(b_) {
            for (int i = 0; i < 4; i++) {
                if (a.getSize(i) && b.getSize(i)) {
                    assert(a.getSize(i) == b.getSize(i),
                           "Can only combine images with matching size\n");
                }
            }
        }

        int getSize(int i) const {
            if (i == 1) {
                return std::max(a.getSize(1), b.getSize(1))*2;
            } else {
                return std::max(a.getSize(i), b.getSize(i));
            }
        }

        float operator()(int x, int y, int t, int c) const {
            if (y & 1) return b(x, y/2, t, c);
            else return a(x, y/2, t, c);
        }

        struct Iter {
            const typename A::Iter a;
            const typename B::Iter b;
            const bool useA;
            Iter() : useA(false) {}
            // Second argument exists to distinguish between the two
            // constructors when A and B have the same type
            Iter(const typename A::Iter &a_, bool) : a(a_), useA(true) {}
            Iter(const typename B::Iter &b_) : b(b_), useA(false) {}
            float operator[](int x) const {
                if (useA) return a[x];
                else return b[x];
            }
            Vec::type vec(int x) const {                
                if (useA) return a.vec(x);
                else return b.vec(x);
            }
        };
 
        Iter scanline(int x, int y, int t, int c, int width) const {
            if (y & 1) return Iter(b.scanline(x, y/2, t, c, width));
            else return Iter(a.scanline(x, y/2, t, c, width), false);
        }

        void prepare(int x, int y, int t, int c, 
                     int width, int height, int frames, int channels) const {
            a.prepare(x, (y+1)/2, t, c, width, height/2, frames, channels);
            b.prepare(x, y/2, t, c, width, height/2, frames, channels);
        }
    };

    template<typename A, typename B>
    _InterleaveY<typename A::Lazy, typename B::Lazy>
    interleaveY(const A &a, const B &b) {
        return _InterleaveY<A, B>(a, b);
    }

    template<typename B>
    _InterleaveY<Const, typename B::Lazy>
    interleaveY(float a, const B &b) {
        return _InterleaveY<Const, B>(Const(a), b);
    }

    template<typename A>
    _InterleaveY<typename A::Lazy, Const>
    interleaveY(const A &a, float b) {
        return _InterleaveY<A, Const>(a, Const(b));
    }

    inline _InterleaveY<Const, Const>
    interleaveY(float a, float b) {
        return _InterleaveY<Const, Const>(Const(a), Const(b));
    }

    template<typename A, typename B>
    struct _InterleaveT {
        typedef _InterleaveT<typename A::Lazy, typename B::Lazy> Lazy;
        typename Handle<A>::type a;
        typename Handle<B>::type b;
        
        _InterleaveT(const A &a_, const B &b_) : a(a_), b(b_) {
            for (int i = 0; i < 4; i++) {
                if (a.getSize(i) && b.getSize(i)) {
                    assert(a.getSize(i) == b.getSize(i),
                           "Can only combine images with matching size\n");
                }
            }
        }

        int getSize(int i) const {
            if (i == 2) {
                return std::max(a.getSize(2), b.getSize(2))*2;
            } else {
                return std::max(a.getSize(i), b.getSize(i));
            }
        }

        float operator()(int x, int y, int t, int c) const {
            if (t & 1) return b(x, y, t/2, c);
            else return a(x, y, t/2, c);
        }

        typedef typename _InterleaveY<A, B>::Iter Iter;

        Iter scanline(int x, int y, int t, int c, int width) const {
            if (t & 1) return Iter(b.scanline(x, y, t/2, c, width));
            else return Iter(a.scanline(x, y, t/2, c, width), false);
        }

        void prepare(int x, int y, int t, int c, 
                     int width, int height, int frames, int channels) const {
            a.prepare(x, y, (t+1)/2, c, width, height, frames/2, channels);
            b.prepare(x, y, t/2, c, width, height, frames/2, channels);
        }
    };

    template<typename A, typename B>
    _InterleaveT<typename A::Lazy, typename B::Lazy>
    interleaveT(const A &a, const B &b) {
        return _InterleaveT<A, B>(a, b);
    }
    template<typename B>
    _InterleaveT<Const, typename B::Lazy>
    interleaveT(float a, const B &b) {
        return _InterleaveT<Const, B>(Const(a), b);
    }

    template<typename A>
    _InterleaveT<typename A::Lazy, Const>
    interleaveT(const A &a, float b) {
        return _InterleaveT<A, Const>(a, Const(b));
    }

    inline _InterleaveT<Const, Const>
    interleaveT(float a, float b) {
        return _InterleaveT<Const, Const>(Const(a), Const(b));
    }

    template<typename A, typename B>
    struct _InterleaveC {
        typedef _InterleaveC<typename A::Lazy, typename B::Lazy> Lazy;
        typename Handle<A>::type a;
        typename Handle<B>::type b;
        
        _InterleaveC(const A &a_, const B &b_) : a(a_), b(b_) {
            for (int i = 0; i < 4; i++) {
                if (a.getSize(i) && b.getSize(i)) {
                    assert(a.getSize(i) == b.getSize(i),
                           "Can only combine images with matching size\n");
                }
            }
        }

        int getSize(int i) const {
            if (i == 3) {
                return std::max(a.getSize(3), b.getSize(3))*2;
            } else {
                return std::max(a.getSize(i), b.getSize(i));
            }
        }

        float operator()(int x, int y, int t, int c) const {
            if (c & 1) return b(x, y, t, c/2);
            else return a(x, y, t, c/2);
        }

        typedef typename _InterleaveY<A, B>::Iter Iter;

        Iter scanline(int x, int y, int t, int c, int width) const {
            if (c & 1) return Iter(b.scanline(x, y, t, c/2, width));
            else return Iter(a.scanline(x, y, t, c/2, width), false);
        }

        void prepare(int x, int y, int t, int c, 
                     int width, int height, int frames, int channels) const {
            a.prepare(x, y, t, (c+1)/2, width, height, frames, channels/2);
            b.prepare(x, y, t, c/2, width, height, frames, channels/2);
        }
    };

    template<typename A, typename B>
    _InterleaveC<typename A::Lazy, typename B::Lazy>
    interleaveC(const A &a, const B &b) {
        return _InterleaveC<A, B>(a, b);
    }

    template<typename B>
    _InterleaveC<Const, typename B::Lazy>
    interleaveC(float a, const B &b) {
        return _InterleaveC<Const, B>(Const(a), b);
    }

    template<typename A>
    _InterleaveC<typename A::Lazy, Const>
    interleaveC(const A &a, float b) {
        return _InterleaveC<A, Const>(a, Const(b));
    }

    inline _InterleaveC<Const, Const>
    interleaveC(float a, float b) {
        return _InterleaveC<Const, Const>(Const(a), Const(b));
    }

    // A trait to check if something is ok to be cast to a lazy expression type, and if so, how.
    template<typename T>
    struct Lazyable {
        typedef typename T::Lazy t;
    };
    
    // Or consts
    template<>
    struct Lazyable<float> {
        typedef Const t;
    };
    
    template<>
    struct Lazyable<int> {
        typedef Const t;
    };

    template<>
    struct Lazyable<double> {
        typedef Const t;
    };



}

// We need to generate a stupid number of operators to overload.
// Traits can help here, but msvc has quirky behaviour with sfinae, so we'll generate them with macros

// First arg is BinaryOp or Cmp
// Second arg is the Symbol (e.g. +)
// Third arg is the Lazy::Vec struct that does the operation (e.g. Add)
// There are three macros - the one that takes two lazy args,
// and the ones where one of the args is a numeric const (float, int, double).
// In this second case, the fourth arg is the type of the numeric const.
#define MAKE_OP_LL(T, S, N)                                     \
    template<typename A, typename B>                            \
    Lazy::T<typename A::Lazy, typename B::Lazy, Lazy::Vec::N>   \
    operator S(const A &a, const B &b) {                        \
        return Lazy::T<A, B, Lazy::Vec::N>(a, b);               \
    }

#define MAKE_OP_CL(T, S, N, CT)                                         \
    template<typename B>                                                \
    Lazy::T<Lazy::Const, typename B::Lazy, Lazy::Vec::N>                \
    operator S(CT a, const B &b) {                                      \
        return Lazy::T<Lazy::Const, B, Lazy::Vec::N>(Lazy::Const(a), b); \
    }

#define MAKE_OP_LC(T, S, N, CT)                                         \
    template<typename A>                                                \
    Lazy::T<typename A::Lazy, Lazy::Const, Lazy::Vec::N>                \
    operator S(const A &a, CT b) {                                      \
        return Lazy::T<A, Lazy::Const, Lazy::Vec::N>(a, Lazy::Const(b)); \
    }

// Make the full set of operator overloads for a given operator
#define MAKE_OP(T, S, N)                        \
    MAKE_OP_LL(T, S, N)                         \
    MAKE_OP_LC(T, S, N, float)                  \
    MAKE_OP_LC(T, S, N, double)                 \
    MAKE_OP_LC(T, S, N, int)                    \
    MAKE_OP_CL(T, S, N, float)                  \
    MAKE_OP_CL(T, S, N, double)                 \
    MAKE_OP_CL(T, S, N, int)                    \
     

MAKE_OP(BinaryOp, +, Add)
MAKE_OP(BinaryOp, -, Sub)
MAKE_OP(BinaryOp, *, Mul)
MAKE_OP(Cmp, >, GT)
MAKE_OP(Cmp, <, LT)
MAKE_OP(Cmp, >=, GE)
MAKE_OP(Cmp, <=, LE)
MAKE_OP(Cmp, ==, EQ)
MAKE_OP(Cmp, !=, NEQ)

// Unary negation is a special case
template<typename A>
Lazy::BinaryOp<Lazy::Const, typename A::Lazy, Lazy::Vec::Sub>
operator-(const A &a) {
    return Lazy::Const(0.0f) - a;
}

// Division by a scalar is another special case
MAKE_OP_LL(BinaryOp, /, Div)
MAKE_OP_CL(BinaryOp, /, Div, float)
MAKE_OP_CL(BinaryOp, /, Div, double)
MAKE_OP_CL(BinaryOp, /, Div, int)

template<typename A>
Lazy::BinaryOp<typename A::Lazy, Lazy::Const, Lazy::Vec::Mul>
operator/(const A &a, float b) {
    return a * (1.0f/b);
}
template<typename A>
Lazy::BinaryOp<typename A::Lazy, Lazy::Const, Lazy::Vec::Mul>
operator/(const A &a, int b) {
    return a * (1.0f/b);
}
template<typename A>
Lazy::BinaryOp<typename A::Lazy, Lazy::Const, Lazy::Vec::Mul>
operator/(const A &a, double b) {
    return a * (1.0f/b);
}

#include "footer.h"

#endif
