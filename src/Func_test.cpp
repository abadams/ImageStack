#include "ImageStack.h"
#include "Func.h"

using namespace ImageStack;
using namespace ImageStack::Lazy;

int main(int argc, char **argv) {
    start();

    try {
        Image in = Load::apply("in.tmp");       
        Image left(in.width, in.height, in.frames, in.channels);
        left.set(shiftX(zeroBoundary(in), 1));
        Image right(in.width, in.height, in.frames, in.channels);
        right.set(shiftX(zeroBoundary(in), -1));
        Image up = interleaveX(3*in + shiftX(zeroBoundary(in), 1), 3*in + shiftX(zeroBoundary(in), -1))/4;
        Save::apply(up, "up.tmp");
        Save::apply(left, "left.tmp");
        Save::apply(right, "right.tmp");
    } catch (Exception &e) {
        printf("Failure: %s\n", e.message);
    }

    return 0;
}





