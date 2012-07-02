#ifndef IMAGESTACK_MACROS_H
#define IMAGESTACK_MACROS_H

#include <stdarg.h>

#include "header.h"

#ifdef WIN32
#undef min
#undef max
#endif

template<typename T>
inline T max(const T &a, const T &b) {
    if (a > b) { return a; }
    return b;
}

template<typename T>
inline T min(const T &a, const T &b) {
    if (a < b) { return a; }
    return b;
}

template<typename T>
inline T max(const T &a, const T &b, const T &c) {
    if (a > b && a > c) { return a; }
    if (b > c) { return b; }
    return c;
}

template<typename T>
inline T min(const T &a, const T &b, const T &c) {
    if (a < b && a < c) { return a; }
    if (b < c) { return b; }
    return c;
}

template<typename T>
inline T clamp(const T &a, const T &b, const T &c) {
    if (a < b) { return b; }
    if (a > c) { return c; }
    return a;
}

// Equal to within 1/100. Useful for testing
static inline bool nearlyEqual(float a, float b) {
    if (fabs(b) < 1) return fabs(a-b) < 0.01;
    return fabs(a/b - 1.0) < 0.01;
}

#ifndef M_PI
#define M_PI 3.14159265
#endif

#ifndef E
#define E 2.7182818284590451
#endif

#define INF (std::numeric_limits<float>::infinity())

static inline float sinc(float x) {
    return x == 0 ? 1 : sinf(M_PI * x) / (M_PI * x);
}

// Map high dynamic range values to [0, 255], so that 0->0, and 1->255
static inline unsigned char HDRtoLDR(float x) {
    if (x < 0) { return 0; }
    if (x > 1) { return 255; }
    return (unsigned char)(x * 255.0f + 0.49999f);
}

// Map low dynamic range values to [0, 1] inclusive.
static inline float LDRtoHDR(unsigned char x) {
    return x * (1.0f/255);
}

static inline float LDR16toHDR(unsigned short x) {
    return x * (1.0f/65535);
}

#include <math.h>

#include "footer.h"
#endif
