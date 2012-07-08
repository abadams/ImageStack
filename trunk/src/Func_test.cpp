#include "ImageStack.h"
#include "Func.h"

using namespace ImageStack;
using namespace ImageStack::Lazy;

int main(int argc, char **argv) {
    start();

    try {
        Image in = Load::apply("in.tmp");       
        Image up(in.width*2, in.height, in.frames, in.channels);

        double t1 = currentTime();
        for (int i = 0; i < 100; i++) {
            up.set(interleaveX(in, in)*0.75f + 
                   interleaveX(shiftX(zeroBoundary(in), 1), 
                               shiftX(zeroBoundary(in), -1))*0.25f);
        }
        double t2 = currentTime();

        //Image up(in.width, in.height, in.frames, in.channels);
        for (int i = 0; i < 100; i++) {
            
            Image inM = in.selectColumns(1, in.width-2);
            Image inL = in.selectColumns(0, in.width-2);
            Image inR = in.selectColumns(2, in.width-2);
            up.selectColumns(2, 2*in.width-4).set(
                interleaveX(inM, inM)*0.75f + interleaveX(inL, inR)*0.25f
                );
        }
        double t3 = currentTime();
        
        for (int i = 0; i < 100; i++) {
            Func zb = zeroBoundary(interleaveX(in, in));
            up.set((shiftX(zb, -1) + 2*zb + shiftX(zb, 1))/4);
        }
    
        double t4 = currentTime();

        printf("%f %f %f\n", t2-t1, t3-t2, t4-t3);
        
        Save::apply(up, "up.tmp");
    } catch (Exception &e) {
        printf("Failure: %s\n", e.message);
    }

    return 0;
}





