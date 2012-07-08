#include "ImageStack.h"
#include "Func.h"

using namespace ImageStack;
using namespace ImageStack::Lazy;

int main(int argc, char **argv) {
    start();

    try {
        Image in = Load::apply("in.tmp");       
        Image up(in.width, in.height, in.frames, in.channels);

        double t1 = currentTime();
        printf("Inline...\n");
        for (int i = 0; i < 10; i++) {
            // blur the square root
            // inline
            auto tmp = zeroBoundary((in));
            auto bx = shiftX(tmp, -1) + 2*tmp + shiftX(tmp, 1);
            auto by = shiftY(bx, -1) + 2*bx + shiftY(bx, 1);
            up.set(by);
        }
        double t2 = currentTime();

        printf("Root...\n");
        for (int i = 0; i < 10; i++) {
            // root
            Image tmp = (in);
            Image bx = shiftX(zeroBoundary(tmp), -1) + 2*tmp + shiftX(zeroBoundary(tmp), 1);
            Image by = shiftY(zeroBoundary(bx), -1) + 2*bx + shiftY(zeroBoundary(bx), 1);
            up = by;
        }
        double t3 = currentTime();

        printf("Chunk...\n");        
        for (int i = 0; i < 10; i++) {
            // chunk
            Func tmp = zeroBoundary((in));
            Func bx = shiftX(tmp, -1) + 2*tmp + shiftX(tmp, 1);
            Func by = shiftY(bx, -1) + 2*bx + shiftY(bx, 1);
            up.set(by);
        }
        double t4 = currentTime();

        printf("%f %f %f\n", t2-t1, t3-t2, t4-t3);
        
        Save::apply(up, "up.tmp");
    } catch (Exception &e) {
        printf("Failure: %s\n", e.message);
    }

    return 0;
}





