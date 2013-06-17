#ifndef IMAGESTACK_CONTROL_H
#define IMAGESTACK_CONTROL_H
namespace ImageStack {

class Loop : public Operation {
public:
    void help();
    bool test() {return true;}
    void parse(vector<string> args);
};

class Pause : public Operation {
public:
    void help();
    bool test() {return true;}
    void parse(vector<string> args);
};

class Time : public Operation {
public:
    void help();
    bool test() {return true;}
    void parse(vector<string> args);
};

}
#endif
