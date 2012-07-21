#ifndef IMAGESTACK_EXPR_SCALAR_H
#define IMAGESTACK_EXPR_SCALAR_H

#include "header.h"

// This file contains scalar definitions of various operators used by Expr.h

namespace Scalar {
    
    // Arithmetic binary operators
    struct Add {
        static float scalar_f(float a, float b) {return a + b;}
        static int scalar_i(int a, int b) {return a + b;}
        
        // Interval arithmetic versions
        template<typename T>
        static std::pair<T, T> interval(std::pair<T, T> a, std::pair<T, T> b) {
            return std::make_pair(a.first+b.first, a.second+b.second);
        }
        
    };
    struct Sub {
        static float scalar_f(float a, float b) {return a - b;}
        static int scalar_i(int a, int b) {return a - b;}
        
        template<typename T>
        static std::pair<T, T> interval(std::pair<T, T> a, std::pair<T, T> b) {
            return std::make_pair(a.first - b.second, a.second - b.first);
        }
    };
    struct Mul {
        static float scalar_f(float a, float b) {return a * b;}
        static int scalar_i(int a, int b) {return a * b;}
        
        template<typename T>
        static std::pair<T, T> interval(std::pair<T, T> a, std::pair<T, T> b) {
            T v1 = a.first * b.first;
            T v2 = a.first * b.second;
            T v3 = a.second * b.first;
            T v4 = a.second * b.second;
            return std::make_pair(
                std::min(std::min(v1, v2), std::min(v3, v4)), 
                std::max(std::max(v1, v2), std::max(v3, v4)));
        }

    };

    struct UnboundedDivisionException {};

    struct Div {
        static float scalar_f(float a, float b) {return a / b;}
        static int scalar_i(int a, int b) {return a / b;}

        template<typename T>
        static std::pair<T, T> interval(std::pair<T, T> a, std::pair<T, T> b) {
            if (b.first <= 0 && b.second >= 0) throw UnboundedDivisionException();
            T v1 = a.first / b.first;
            T v2 = a.first / b.second;
            T v3 = a.second / b.first;
            T v4 = a.second / b.second;
            return std::make_pair(
                std::min(std::min(v1, v2), std::min(v3, v4)), 
                std::max(std::max(v1, v2), std::max(v3, v4)));
        }

    };
    struct Min {
        static float scalar_f(float a, float b) {return std::min(a, b);}
        static int scalar_i(int a, int b) {return std::min(a, b);}

        template<typename T>
        static std::pair<T, T> interval(std::pair<T, T> a, std::pair<T, T> b) {
            return std::make_pair(std::min(a.first, b.first),
                                  std::min(a.second, b.second));
        }
    };
    struct Max {
        static float scalar_f(float a, float b) {return std::max(a, b);}
        static int scalar_i(int a, int b) {return std::max(a, b);}

        template<typename T>
        static std::pair<T, T> interval(std::pair<T, T> a, std::pair<T, T> b) {
            return std::make_pair(std::max(a.first, b.first),
                                  std::max(a.second, b.second));
        }
    };

    // Comparisons
    struct GT {
        static bool scalar_f(float a, float b) {return a > b;}
        static bool scalar_i(int a, int b) {return a > b;}
    };
    struct LT {
        static bool scalar_f(float a, float b) {return a < b;}
        static bool scalar_i(int a, int b) {return a < b;}
    };
    struct GE {
        static bool scalar_f(float a, float b) {return a >= b;}
        static bool scalar_i(int a, int b) {return a >= b;}
    };
    struct LE {
        static bool scalar_f(float a, float b) {return a <= b;}
        static bool scalar_i(int a, int b) {return a <= b;}
    };
    struct EQ {
        static bool scalar_f(float a, float b) {return a == b;}
        static bool scalar_i(int a, int b) {return a == b;}
    };
    struct NEQ {
        static bool scalar_f(float a, float b) {return a != b;}
        static bool scalar_i(int a, int b) {return a != b;}
    };        

    // Unary ops
    struct Ceil {
        static float scalar_f(float a) {return ceilf(a);}

        static std::pair<float, float> interval(std::pair<float, float> a) {
            return make_pair(scalar_f(a.first), scalar_f(a.second));
        }
    };

    struct Floor {
        static float scalar_f(float a) {return floorf(a);}

        static std::pair<float, float> interval(std::pair<float, float> a) {
            return make_pair(scalar_f(a.first), scalar_f(a.second));
        }
    };

    struct Sqrt {
        static float scalar_f(float a) {return sqrtf(a);}

        static std::pair<float, float> interval(std::pair<float, float> a) {
            return make_pair(scalar_f(a.first), scalar_f(a.second));
        }
    };
}

#include "footer.h"

#endif
