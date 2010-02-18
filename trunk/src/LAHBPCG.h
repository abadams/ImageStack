#ifndef LAHBPCG_H
#define LAHBPCG_H

class LAHBPCG : public Operation {
  public:
    void help();
    void parse(vector<string> args);
  
    static Image apply(Window d, Window gx, Window gy, Window w, Window sx, Window sy, int max_iter, float tol);
  private:
};

#endif
