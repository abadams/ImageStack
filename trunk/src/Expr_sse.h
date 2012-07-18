#ifndef IMAGESTACK_EXPR_SSE_H
#define IMAGESTACK_EXPR_SSE_H

#include <immintrin.h>

#include "header.h"

namespace Vec {
    typedef __m128 type;
    const int width = 4;
    
    inline type broadcast(float v) {
        return _mm_set1_ps(v);
    }
    
    inline type set(float a, float b, float c, float d, float e = 0, float f = 0, float g = 0, float h = 0) {
        return _mm_set_ps(d, c, b, a);
    }
    
    inline type zero() {
        return _mm_setzero_ps();
    }
    
    // Arithmetic binary operators
    struct Add : public ImageStack::Scalar::Add {
        static type vec(type a, type b) {return _mm_add_ps(a, b);}
    };
    struct Sub : public ImageStack::Scalar::Sub {
        static type vec(type a, type b) {return _mm_sub_ps(a, b);}
    };
    struct Mul : public ImageStack::Scalar::Mul {
        static type vec(type a, type b) {return _mm_mul_ps(a, b);}
    };
    struct Div : public ImageStack::Scalar::Div {
        static type vec(type a, type b) {return _mm_div_ps(a, b);}
    };
    struct Min : public ImageStack::Scalar::Min {
        static type vec(type a, type b) {return _mm_min_ps(a, b);}
    };
    struct Max : public ImageStack::Scalar::Max {
        static type vec(type a, type b) {return _mm_max_ps(a, b);}
    };
    
    // Comparisons
    struct GT : public ImageStack::Scalar::GT {
        static type vec(type a, type b) {return _mm_cmpgt_ps(a, b);}
    };
    struct LT : public ImageStack::Scalar::LT {
        static type vec(type a, type b) {return _mm_cmplt_ps(a, b);}
    };
    struct GE : public ImageStack::Scalar::GE {
        static type vec(type a, type b) {return _mm_cmpge_ps(a, b);}
    };
    struct LE : public ImageStack::Scalar::LE {
        static type vec(type a, type b) {return _mm_cmple_ps(a, b);}
    };
    struct EQ : public ImageStack::Scalar::EQ {
        static type vec(type a, type b) {return _mm_cmpeq_ps(a, b);}
    };
    struct NEQ : public ImageStack::Scalar::NEQ {
        static type vec(type a, type b) {return _mm_cmpneq_ps(a, b);}
    };
    
    inline type interleave(type a, type b) {
        return _mm_unpacklo_ps(a, b);
    }
    
    inline type subsample(type a, type b) {
        return _mm_shuffle_ps(a, b, (0 << 0) | (2 << 2) | (1 << 4) | (3 << 6));
    }
    
    inline type reverse(type a) {
        return _mm_shuffle_ps(a, a, (3 << 0) | (2 << 2) | (1 << 4) | (0 << 6));
    }
    
#ifdef __SSE4_1__
    // Logical ops
    inline type blend(type a, type b, type mask) {
        return _mm_blendv_ps(a, b, mask);
    }
    
    // Unary ops
    struct Floor : public ImageStack::Scalar::Floor {
        static type vec(type a) {return _mm_floor_ps(a);}
    };
    struct Ceil : public ImageStack::Scalar::Ceil {
        static type vec(type a) {return _mm_ceil_ps(a);}
    };
    struct Sqrt : public ImageStack::Scalar::Sqrt {
        static type vec(type a) {return _mm_sqrt_ps(a);}
    };
#else
    
    inline type blend(type a, type b, type mask) {
        return _mm_or_ps(_mm_and_ps(mask, b),
                         _mm_andnot_ps(mask, a));
    }
    
    struct Floor : public ImageStack::Scalar::Floor {
        static type vec(type a) {
            union {
                float f[width];
                type v;
            } v;
            v.v = a;
            return set(scalar_f(v.f[0]), scalar_f(v.f[1]), scalar_f(v.f[2]), scalar_f(v.f[3]));
        }
    };

    struct Ceil : public ImageStack::Scalar::Ceil {
        static type vec(type a) {
            union {
                float f[width];
                type v;
            } v;
            v.v = a;
            return set(scalar_f(v.f[0]), scalar_f(v.f[1]), scalar_f(v.f[2]), scalar_f(v.f[3]));
        }
    };

    struct Sqrt : public ImageStack::Scalar::Sqrt {
        static type vec(type a) {
            union {
                float f[width];
                type v;
            } v;
            v.v = a;
            return set(scalar_f(v.f[0]), scalar_f(v.f[1]), scalar_f(v.f[2]), scalar_f(v.f[3]));
        }
    };
    
#endif
    
    // Loads and stores
    inline type load(const float *f) {
        return _mm_loadu_ps(f);
    }
    
    inline void store(type a, float *f) {
        _mm_storeu_ps(f, a);
    }
}

#include "footer.h"
#endif
