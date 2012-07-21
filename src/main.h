#ifndef IMAGESTACK_MAIN_H
#define IMAGESTACK_MAIN_H

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <limits>
#include <list>
#include <sstream>
#include <memory>

using ::std::shared_ptr;
using ::std::string;
using ::std::vector;
using ::std::pair;
using ::std::make_pair;
using ::std::map;
using ::std::list;
using ::std::swap;

#ifdef _MSC_VER
#include <windows.h>
#include <float.h>
#define popen _popen
#define pclose _pclose
#define isnan _isnan
#define isfinite _finite
inline float isinf(float x) {
    return (!_finite(x) && !_isnan(x));
}
#else
using ::std::isfinite;
#endif

#ifdef __CYGWIN__
using ::std::isinf;
using ::std::isnan;
#endif

// Some core files that everyone should include
#include "macros.h"
#include "Exception.h"
#include "Operation.h"

namespace ImageStack {
// time since program start in seconds (if using from a library, time since ImageStack::begin)
float currentTime();
}

#include "Image.h"
#include "header.h"

// Below are the data structures and functions available to operations:

// Deal with the stack of images that gives this program its name
Image &stack(size_t index);
void push(Image);
void pop();
void dup();
void pull(size_t);

// Parse ints, floats, chars, and ImageStack commands
int readInt(string);
float readFloat(string);
char readChar(string);
void parseCommands(vector<string>);

// Fire up and shut down imagestack. This populates the operation map,
// and sets a starting time for timing ops.
void start();
void end();

// Generate a uniform random integer within [min, max]
int randomInt(int min, int max);

// Generate a uniform random float within [min, max]
float randomFloat(float min, float max);

// pretty-print some help text, by word wrapping at 79 chars
void pprintf(const char *str);

// Check if two images are nearly equal for unit testing
bool nearlyEqual(Image a, Image b);

// The map of operations, which converts strings to operation
// objects. Only meta-operations like help should need to access this.
extern map<string, Operation *> operationMap;
typedef map<string, Operation *>::iterator OperationMapIterator;

#include "footer.h"
#endif
