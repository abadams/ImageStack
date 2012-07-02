#ifndef IMAGESTACK_STACK_H
#define IMAGESTACK_STACK_H
#include <map>
#include "header.h"

// These operations apply only to the stack, so they have no apply
// method. They also have no unit tests.

class Pop : public Operation {
public:
    void help();
    bool test() {return true;}
    void parse(vector<string> args);
};

class Push : public Operation {
public:
    void help();
    bool test() {return true;}
    void parse(vector<string> args);
};

class Pull : public Operation {
public:
    void help();
    bool test() {return true;}
    void parse(vector<string> args);
};

class Dup : public Operation {
public:
    void help();
    bool test() {return true;}
    void parse(vector<string> args);
};

class Stash : public Operation {
public:
    void help();
    bool test() {return true;}
    void parse(vector<string> args);
    static map<string, Image> stash;
};

#include "footer.h"
#endif
