#ifndef IMAGESTACK_EXCEPTION_H
#define IMAGESTACK_EXCEPTION_H

#include <stdarg.h>

namespace ImageStack {

#define EXCEPTION_LENGTH 1024

class Exception {
public:
    Exception(const char *fmt, va_list arglist);
    Exception(const char *fmt, ...);

    char message[EXCEPTION_LENGTH];
};

void panic(const char *fmt, ...) throw(Exception);
#undef assert
void assert(bool cond, const char *fmt, ...) throw(Exception);

}

#endif
