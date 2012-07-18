#ifndef IMAGESTACK_EXPR_AVX_H
#define IMAGESTACK_EXPR_AVX_H

#include <immintrin.h>

#include "header.h"

namespace Vec {
    typedef __m256 type;
    const int width = 8;
    
    inline type broadcast(float v) {
        return _mm256_set1_ps(v);
    }
    
    inline type set(float a, float b, float c, float d = 0, float e = 0, float f = 0, float g = 0, float h = 0) {
        return _mm256_set_ps(h, g, f, e, d, c, b, a);
    }
    
    inline type zero() {
        return _mm256_setzero_ps();
    }
    
    // Arithmetic binary operators
    struct Add : public ImageStack::Scalar::Add {
        static type vec(type a, type b) {return _mm256_add_ps(a, b);}
    };
    struct Sub : public ImageStack::Scalar::Sub {
        static type vec(type a, type b) {return _mm256_sub_ps(a, b);}
    };
    struct Mul : public ImageStack::Scalar::Mul {
        static type vec(type a, type b) {return _mm256_mul_ps(a, b);}
    };
    struct Div : public ImageStack::Scalar::Div {
        static type vec(type a, type b) {return _mm256_div_ps(a, b);}
    };
    struct Min : public ImageStack::Scalar::Min {
        static type vec(type a, type b) {return _mm256_min_ps(a, b);}
    };
    struct Max : public ImageStack::Scalar::Max {
        static type vec(type a, type b) {return _mm256_max_ps(a, b);}
    };

    // Comparisons
    struct GT : public ImageStack::Scalar::GT {
        static type vec(type a, type b) {return _mm256_cmp_ps(a, b, _CMP_GT_OQ);}
    };
    struct LT : public ImageStack::Scalar::LT {
        static type vec(type a, type b) {return _mm256_cmp_ps(a, b, _CMP_LT_OQ);}
    };
    struct GE : public ImageStack::Scalar::GE {
        static type vec(type a, type b) {return _mm256_cmp_ps(a, b, _CMP_GE_OQ);}
    };
    struct LE : public ImageStack::Scalar::LE {
        static type vec(type a, type b) {return _mm256_cmp_ps(a, b, _CMP_LE_OQ);}
    };
    struct EQ : public ImageStack::Scalar::EQ {
        static type vec(type a, type b) {return _mm256_cmp_ps(a, b, _CMP_EQ_OQ);}
    };
    struct NEQ : public ImageStack::Scalar::NEQ {
        static type vec(type a, type b) {return _mm256_cmp_ps(a, b, _CMP_NEQ_OQ);}
    };

    // Logical ops
    inline type blend(type a, type b, type mask) {
        return _mm256_blendv_ps(a, b, mask);
    }

    inline type interleave(type a, type b) {
        // Given vectors a and b, return a[0] b[0] a[1] b[1] a[2] b[2] a[3] b[3]
        __m256 r_lo = _mm256_unpacklo_ps(a, b);
        __m256 r_hi = _mm256_unpackhi_ps(a, b);
        return _mm256_permute2f128_ps(r_lo, r_hi, 2<<4);
    }

    inline type subsample(type a, type b) {
        // Given vectors a and b, return a[0], a[2], a[4], a[6], b[1], b[3], b[5], b[7]
        type bodd = _mm256_shuffle_ps(b, b, (1 << 0) | (1 << 2) | (3 << 4) | (3 << 6));
        // bodd = b[1] b[1] b[3] b[3] b[5] b[5] b[7] b[7]
        type lo = _mm256_permute2f128_ps(a, bodd, (0 << 0) | (2 << 4));
        // lo = a[0] a[1] a[2] a[3] b[1] b[1] b[3] b[3]
        type hi = _mm256_permute2f128_ps(a, bodd, (1 << 0) | (3 << 4));
        // hi = a[4] a[5] a[6] a[7] b[5] b[5] b[7] b[7]
        type result = _mm256_shuffle_ps(lo, hi, (2 << 2) | (2 << 6));
        // result = a[0] a[2] a[4] a[6] b[1] b[3] b[5] b[7]
        return result;

    }
        
    inline type reverse(type a) {
        // reverse each half, the reverse the halves
        type b = _mm256_shuffle_ps(a, a, (3 << 0) | (2 << 2) | (1 << 4) | (0 << 6));
        return _mm256_permute2f128_ps(b, b, 1);
    }


    // Unary ops
    struct Floor : public ImageStack::Scalar::Floor {
        static type vec(type a) {return _mm256_floor_ps(a);}
    };
    struct Ceil : public ImageStack::Scalar::Ceil {
        static type vec(type a) {return _mm256_ceil_ps(a);}
    };
    struct Sqrt : public ImageStack::Scalar::Sqrt {
        static type vec(type a) {return _mm256_sqrt_ps(a);}
    };

    // Loads and stores
    inline type load(const float *f) {
        return _mm256_loadu_ps(f);
    }

    inline void store(type a, float *f) {
        _mm256_storeu_ps(f, a);
    }
}

#include "footer.h"

#endif
