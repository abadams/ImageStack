#ifndef IMAGESTACK_EXPR_SCALAR_FALLBACK_H
#define IMAGESTACK_EXPR_SCALAR_FALLBACK_H

#include "header.h"

// This file gets included if no sse or avx is available. It provides
// for a vector type that's just a float

namespace Vec {
    typedef float type;
    const int width = 1;

    inline type broadcast(float v) {
        return v;
    }

    inline type zero() {
        return 0;
    }

    // Arithmetic binary operators
    struct Add : public ImageStack::Scalar::Add {
        static type vec(type a, type b) {return scalar_f(a, b);}
    };
    struct Sub : public ImageStack::Scalar::Sub {
        static type vec(type a, type b) {return scalar_f(a, b);}
    };
    struct Mul : public ImageStack::Scalar::Mul {
        static type vec(type a, type b) {return scalar_f(a, b);}
    };
    struct Div : public ImageStack::Scalar::Div {
        static type vec(type a, type b) {return scalar_f(a, b);}
    };
    struct Min : public ImageStack::Scalar::Min {
        static type vec(type a, type b) {return scalar_f(a, b);}
    };
    struct Max : public ImageStack::Scalar::Max {
        static type vec(type a, type b) {return scalar_f(a, b);}
    };

    // Comparisons
    struct GT : public ImageStack::Scalar::GT {
        static bool vec(type a, type b) {return scalar_f(a, b);}
    };
    struct LT : public ImageStack::Scalar::Lt {
        static bool vec(type a, type b) {return scalar_f(a, b);}
    };
    struct GE : public ImageStack::Scalar::GE {
        static bool vec(type a, type b) {return scalar_f(a, b);}
    };
    struct LE : public ImageStack::Scalar::LE {
        static bool vec(type a, type b) {return scalar_f(a, b);}
    };
    struct EQ : public ImageStack::Scalar::EQ {
        static bool vec(type a, type b) {return scalar_f(a, b);}
    };
    struct NEQ : public ImageStack::Scalar::NEQ {
        static bool vec(type a, type b) {return scalar_f(a, b);}
    };

    // Logical ops
    inline type blend(bool mask, type a, type b) {
        return (mask ? b : a);
    }

    inline type interleave(type a, type b) {
        return a;
    }

    inline type subsample(type a, type b) {
        return a;
    }

    inline type reverse(type a) {
        return a;
    }

    // Unary ops
    struct Floor : public ImageStack::Scalar::Floor {
        static type vec(type a) {return scalar_f(a);}
    };
    struct Ceil : public ImageStack::Scalar::Ceil {
        static type vec(type a) {return scalar_f(a);}
    };
    struct Sqrt : public ImageStack::Scalar::Sqrt {
        static type vec(type a) {return scalar_f(a);}
    };

    // Loads and stores
    inline type load(const float *f) {
        return *f;
    }

    inline void store(type a, float *f) {
        *f = a;
    }
}


#include "footer.h"

#endif
