// An example plugin to be used with the -plugin operator.
//
// On OS X you compile it like this:
// g++ -c plugin.cpp -I ../src
// ld -dylib plugin.o -o plugin.so -undefined dynamic_lookup
//
// On linux you compile it like this:
// g++ -shared plugin.cpp -o plugin.so -I ../src -fPIC
//
// And use it like this:
// ImageStack -load ../pics/dog1.jpg -plugin ./plugin.so -help foo -foo -display

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

    // Foo doubles an image and then squares it
    static Image apply(Window im) {
        Image out(im);

        // Call an existing imagestack operator
        Add::apply(out, out);

        // Do some custom stuff
        for (int t = 0; t < out.frames; t++) {
            for (int y = 0; y < out.height; y++) {
                for (int x = 0; x < out.width; x++) {
                    for (int c = 0; c < out.channels; c++) {
                        out(x, y, t)[c] *= out(x, y, t)[c];
                    }
                }
            }
        }

        return out;
    }

    void help() {
        pprintf("Foo help\n");
    }
};

extern "C" void init_imagestack_plugin(map<string, Operation *> &operationMap) {
    printf("Initializing ImageStack plugin\n");
    operationMap["-foo"] = new Foo();
    printf("Added foo to the operation map\n");
}
