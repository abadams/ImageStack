#ifndef IMAGESTACK_PROJECTION_H
#define IMAGESTACK_PROJECTION_H

class Sinugram : public Operation {
  public:
    void help();
    void parse(vector<string> args);
    static Image apply(Window im, int directions);
};

#endif
