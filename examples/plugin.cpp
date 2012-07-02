// An example plugin to be used with the -plugin operator.
//
// On OS X you compile it like this:
// g++ -c plugin.cpp -I ../src
// ld -dylib plugin.o -o plugin.so -undefined dynamic_lookup
//
// On linux you compile it like this:
// g++ -std=gnu++0x -shared plugin.cpp -o plugin.so -I ../src -fPIC
//
// And use it like this:
// ImageStack -load ../pics/dog1.jpg -plugin ./plugin.so -help foo -test foo -foo -display

#include <ImageStack.h>

using namespace ImageStack;

// here's our new operation to add
class Foo : public Operation {
public:
    void parse(vector<string> args) {
        printf("Called foo with %zu args\n", args.size());
        Image im = apply(stack(0));
        pop();
        push(im);
    }

    // Foo doubles an image then squares it, then flips it
    static Image apply(Image im) {
        Image out(im);

        // Use an image expression
        out += out;

        // Do some custom stuff
        for (int t = 0; t < out.frames; t++) {
            for (int y = 0; y < out.height; y++) {
                for (int x = 0; x < out.width; x++) {
                    for (int c = 0; c < out.channels; c++) {
                        out(x, y, t, c) *= out(x, y, t, c);
                    }
                }
            }
        }

        // Use an imagestack operator
        Flip::apply(out, 'x');

        return out;
    }

    void help() {
        pprintf("-foo doubles an image, squares it, and flips it.\n");
    }

    bool test() {
        // Apply foo to a linear ramp
        Image ramp(128, 128, 1, 1);
        Lazy::X x; Lazy::Y y;
        ramp.set(x + 3*y);
        Image out = Foo::apply(ramp);

        // For a ramp the output should be...
        Image correct(128, 128, 1, 1);
        correct.set((ramp.width-1-x) + 3*y);
        correct *= 2;
        correct *= correct;
        
        return nearlyEqual(out, correct);
    }
};


extern "C" void init_imagestack_plugin(map<string, Operation *> &operationMap) {
    printf("Initializing ImageStack plugin\n");
    operationMap["-foo"] = new Foo();
    printf("Added foo to the operation map\n");
}
