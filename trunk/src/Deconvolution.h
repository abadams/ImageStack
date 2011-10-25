#ifndef DECONVOLUTION_H
#define DECONVOLUTION_H
#include "header.h"

/* I am aware that Deconvolve operation exists in Convolve.h, but I did not want to
   touch those, so I am creating a separate file and operation. */
class Deconvolution : public Operation {
 public:
  void help();
  void parse(vector<string> args);
  static Image applyCho2009(Window im, Window kernel);
  static Image applyShan2008(Window im, Window kernel);
 private:
  static Image applyPadding(Window im);
};

#include "footer.h"
#endif
