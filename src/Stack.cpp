#include "main.h"
#include "Stack.h"
#include <list>
#include <map>
#include "header.h"

void Pop::help() {
    pprintf("-pop removes the top image from the stack\n\n"
            "Usage: ImageStack -load a.tga -load b.tga -save b.gif -pop -save a.gif\n");
}

void Pop::parse(vector<string> args) {
    assert(args.size() == 0, "-pop takes no arguments\n");
    pop();
}

void Push::help() {
    pprintf("-push adds a new zeroed image to the top of the stack. With no"
            " arguments it matches the dimensions of the current image. With 4"
            " arguments (width, height, frames, and channels) it creates an image of"
            " that size. Given three arguments frames defaults to 1, and the"
            " arguments are taken as width, height, and channels.\n"
            "\n"
            "Usage: ImageStack -load a.tga -push -add -scale 0.5 -multiply -save out.tga\n"
            "       ImageStack -push 1024 1024 1 3 -offset 0.5 -save gray.tga\n\n");
}

void Push::parse(vector<string> args) {
    if (args.size() == 0) {
        push(Image(stack(0).width, stack(0).height, stack(0).frames, stack(0).channels));
    } else if (args.size() == 3) {
        push(Image(readInt(args[0]), readInt(args[1]), 1, readInt(args[2])));
    } else if (args.size() == 4) {
        push(Image(readInt(args[0]), readInt(args[1]), readInt(args[2]), readInt(args[3])));
    } else {
        panic("-push takes zero, three, or four arguments\n");
    }
}

void Pull::help() {
    pprintf("-pull brings a buried stack element to the top. -pull 0 does nothing. -pull 1"
            " brings up the second stack element, and so on. -pull can also be given"
            " named arguments to retrieve images that have been stashed with"
            " -stash\n"
            "\n"
            "Usage: ImageStack -load a.tga -load b.tga -save b.gif -pull 1 -save a.gif\n\n");
}

void Pull::parse(vector<string> args) {
    assert(args.size() == 1, "-pull takes 1 argument\n");

    if ('1' <= args[0][0] && args[0][0] <= '9') {
        int depth = readInt(args[0]);
        assert(depth > 0, "-pull only makes sense on strictly positive depths\n");
        pull(depth);
    } else {
        map<string, Image>::iterator iter = Stash::stash.find(args[0]);
        assert(iter != Stash::stash.end(),
               "Image with name %s was not found in the stash\n",
               args[0].c_str());
        Image im = iter->second;
        push(im);
        Stash::stash.erase(iter);
    }
}

void Dup::help() {
    pprintf("-dup duplicates an image and pushes it on the stack. Given no argument"
            " it duplicates the top image in the stack. Given a numeric argument it"
            " duplicates that image from down the stack. Given a string argument it"
            " duplicates an image that was stashed with -stash using that name\n"
            "\n"
            "Usage: ImageStack -load a.tga -dup -scale 0.5 -save a_small.tga\n"
            "                  -pop -scale 2 -save a_big.tga\n\n");

}

void Dup::parse(vector<string> args) {
    if (args.size() == 0) {
        dup();
    } else {
        assert(args.size() == 1, "-dup takes zero or one arguments\n");
        if ('0' <= args[0][0] && args[0][0] <= '9') {
            int depth = readInt(args[0]);
            push(stack(depth).copy());
        } else {
            map<string, Image>::iterator iter = Stash::stash.find(args[0]);
            assert(iter != Stash::stash.end(),
                   "Image with name %s was not found in the stash\n",
                   args[0].c_str());
            Image im = iter->second;
            push(im.copy());
        }
    }
}

map<string, Image> Stash::stash;

void Stash::help() {
    pprintf("-stash removes the top image from the stack and gives it a name. It"
            " can be retrieved using -dup or -pull using its name as the argument\n"
            "\n"
            "Usage: ImageStack ... -stash foo ... -dup foo ... -pull foo\n");
}

void Stash::parse(vector<string> args) {
    assert(args.size() == 1, "-stash takes one argument\n");
    Image im = stack(0);
    pop();
    stash[args[0]] = im;
}

#include "footer.h"
