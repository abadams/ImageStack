#include "ImageStack.h"
#include "Func.h"

using namespace ImageStack;
using namespace ImageStack::Lazy;

int main(int argc, char **argv) {
    Image foo(100, 100, 1, 1);
    Image bar(100, 100, 1, 1);
    Image small(30, 30, 1, 1);

    start();

    try {
        // warm up
        Noise::apply(small, 0, 1);
        foo.set(zeroBoundary(sqrt(small)));
        foo.set(shift(zeroBoundary(sqrt(small)), 5, 5));
        Save::apply(foo, "foo.tmp");

        Image dog = Load::apply("pics/dog1.jpg");
        
        printf("Interleaving foo and dog...\n");
        asm("# begin");
        //bar = interleaveC(interleaveT(interleaveY(interleaveX(dog*1, dog*2), interleaveX(dog*4, dog*8)), 0.5f), 0.0f);
        bar = interleaveC(dog, 0.0f);
        asm("# end");
        printf("Done\n");
        Save::apply(bar, "bar.tmp");

        for (int i = 0; i < 10; i++) {
            Image dx(dog.width, dog.height, dog.frames, dog.channels);
            double t1 = currentTime();
            for (int j = 0; j < 10; j++) {
                auto dz = zeroBoundary(dog);
                dx.set(dog + shiftX(dz, -1) + shiftX(dz, 1) + shiftY(dz, 1) + shiftY(dz, -1));
            }
            int w = dog.width, h = dog.height;
            double t2 = currentTime();
            for (int j = 0; j < 10; j++) {
              dx.region(1, 1, 0, 0, w-2, h-2, 1, 3).set(
                  dog.region(1, 1, 0, 0, w-2, h-2, 1, 3) + 
                  dog.region(1, 2, 0, 0, w-2, h-2, 1, 3) + 
                  dog.region(1, 0, 0, 0, w-2, h-2, 1, 3) + 
                  dog.region(2, 1, 0, 0, w-2, h-2, 1, 3) + 
                  dog.region(0, 1, 0, 0, w-2, h-2, 1, 3));
            }
            double t3 = currentTime();
            printf("%f %f\n", (t2-t1)*1000, (t3-t2)*1000);
        }            

    /*printf("Defining g\n");
        Func g = 3*foo;
        printf("Defining f\n");
        Func f = g - 17;
        printf("Setting foo to f * 81\n");
        bar.set(f * 81 + Shift<ZeroBoundary<Func> >(ZeroBoundary<Func>(g), -1, -1, 0, 0));
        printf("Done\n");
    */

    } catch (Exception &e) {
        printf("Failure: %s\n", e.message);
    }

    return 0;
}





