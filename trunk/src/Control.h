#ifndef IMAGESTACK_CONTROL_H
#define IMAGESTACK_CONTROL_H


class Loop : public Operation {
  public:
    void help();
    void parse(vector<string> args);
};

class Pause : public Operation {
  public:
    void help();
    void parse(vector<string> args);
};

class Time : public Operation {
  public:
    void help();
    void parse(vector<string> args);
};

#endif
