#ifndef IMAGESTACK_PREDICTION_H
#define IMAGESTACK_PREDICTION_H

class Inpaint : public Operation {
  public:
    void help();
    void parse(vector<string> args);
    static Image apply(Window im, Window mask);
};

#endif
