#include "main.h"
#include "Plugin.h"
#include "header.h"

#ifdef __WIN32__
#else
#include "dlfcn.h"
#endif

typedef void (*init_imagestack_plugin)(map<string, Operation *> &);

#ifdef __WIN32__

void Plugin::parse(vector<string> args) {
    HMODULE handle = LoadLibrary(args[0].c_str());
    assert(handle, "Could not open %s as a shared library\n", args[0].c_str());

    FARPROC function = GetProcAddress(handle, "init_imagestack_plugin");
    assert(function, "Could not find symbol init_imagestack_plugin (with C linkage) in %s\n", args[0].c_str());

    // Pass control to the shared library to inject crap into the operation map
    ((init_imagestack_plugin)function)(operationMap);
}

#else

void Plugin::parse(vector<string> args) {
    void *handle = dlopen(args[0].c_str(), RTLD_LAZY);
    assert(handle, "Could not open %s as a shared library: %s\n", args[0].c_str(), dlerror());

    void *function = dlsym(handle, "init_imagestack_plugin");
    assert(function, "Could not find symbol init_imagestack_plugin (with C linkage) in %s\n", args[0].c_str());

    // Pass control to the shared library to inject crap into the operation map
    ((init_imagestack_plugin)function)(operationMap);
}

#endif

bool Plugin::test() {
    // The build dependencies to test this are too much.
    return true;
}

void Plugin::help() {
    pprintf("-plugin loads a shared object that can add new operations to"
            " ImageStack. It does this by calling the function"
            " init_imagestack_plugin, which should exist in the shared object with"
            " C linkage. To this function it passes the ImageStack operation map,"
            " into which the plugin may inject new operations.\n\n"
            "Usage: ImageStack -plugin foo.so -foo 1 2 3\n");
}

#include "footer.h"
