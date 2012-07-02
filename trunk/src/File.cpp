#include "main.h"
#include "File.h"
#include "Color.h"
#include "Stack.h"
#include "Arithmetic.h"
#include "Statistics.h"
#include "Filter.h"
#include "header.h"

namespace {
// used for picking file formats
bool suffixMatch(string filename, string suffix) {
    if (suffix.size() > filename.size()) { return false; }
    int offset = (int)filename.size() - (int)suffix.size();
    for (size_t i = 0; i < suffix.size(); i++) {
        if (tolower(filename[i + offset]) != tolower(suffix[i])) { return false; }
    }
    return true;
}


struct TempFile {
    string name;
    TempFile() : name("_test") {}
    TempFile(string name_) : name(name_) {}
    ~TempFile() {
        remove(name.c_str());
    }
};

// Used to help testing. Saves and loads and checks the result is what you saved.
bool testFormat(Image im, string fmt) {
    printf("%s ", fmt.c_str());
    fflush(stdout);
    TempFile t(string("_test") + "." + fmt);
    Save::apply(im, t.name);
    Image b = Load::apply(t.name);
    return nearlyEqual(im, b);
}


};

void Load::help() {
    pprintf("-load loads a file and places it on the top of the stack. ImageStack"
            " can load the following file formats:\n");

    printf("\n");
    FileTMP::help();
    printf("\n");
    FileHDR::help();
    printf("\n");
    FileJPG::help();
    printf("\n");
    FileEXR::help();
    printf("\n");
    FilePNG::help();
    printf("\n");
    FilePPM::help();
    printf("\n");
    FileTGA::help();
    printf("\n");
    FileTIFF::help();
    printf("\n");
    FileWAV::help();
    printf("\n");
    FileFLO::help();
    printf("\n");
    FileCSV::help();
    printf("\n");
    FilePBA::help();
    printf("\n");
    printf("Usage: ImageStack -load foo.jpg\n\n");
}

bool Load::test() {
    // let's take a file around the block and see how it survives
    Image a(123, 234, 3, 3);
    Noise::apply(a, 0, 1);
    FastBlur::apply(a, 1, 1, 1);
    Quantize::apply(a, 1.0/256);

    testFormat(a, "tmp");

    // tmp is the only multi-frame format, so now we switch to a single frame
    a = a.frame(0);

#ifndef NO_JPG
    if (!testFormat(a, "jpg")) return false;
#endif
#ifndef NO_PNG
    if (!testFormat(a, "png")) return false;
#endif
#ifndef NO_TIFF
    if (!testFormat(a, "tiff")) return false;
#endif
#ifndef NO_OPENEXR
    if (!testFormat(a*5, "exr")) return false;
#endif

    if (!testFormat(a, "tga")) return false;
    if (!testFormat(a*5, "hdr")) return false;
    if (!testFormat(a, "ppm")) return false;

    // Now a two-channel format
    a = a.selectChannels(0, 2);
    if (!testFormat(a, "flo")) return false;

    // A two-channel one-row format
    if (!testFormat(a.row(0), "wav")) return false;

    // Now some grayscale formats
    a = a.channel(0);
    if (!testFormat(a, "csv")) return false;
    if (!testFormat(a, "pba")) return false;
    if (!testFormat(a, "pgm")) return false;
    printf("\n");

    return true;
}

void Load::parse(vector<string> args) {
    assert(args.size() == 1, "-load takes exactly 1 argument\n");
    push(apply(args[0]));
}


Image Load::apply(string filename) {
    if (suffixMatch(filename, ".tmp")) {
        return FileTMP::load(filename);
    } else if (suffixMatch(filename, ".hdr")) {
        return FileHDR::load(filename);
    } else if (suffixMatch(filename, ".jpg") ||
               suffixMatch(filename, ".jpeg")) {
        return FileJPG::load(filename);
    } else if (suffixMatch(filename, ".exr")) {
        return FileEXR::load(filename);
    } else if (suffixMatch(filename, ".png")) {
        return FilePNG::load(filename);
    } else if (suffixMatch(filename, ".tga")) {
        return FileTGA::load(filename);
    } else if (suffixMatch(filename, ".wav")) {
        return FileWAV::load(filename);
    } else if (suffixMatch(filename, ".ppm") ||
               suffixMatch(filename, ".pgm")) {
        return FilePPM::load(filename);
    } else if (suffixMatch(filename, ".tiff") ||
               suffixMatch(filename, ".tif")) {
        return FileTIFF::load(filename);
    } else if (suffixMatch(filename, ".flo")) {
        return FileFLO::load(filename);
    } else if (suffixMatch(filename, ".csv")) {
        return FileCSV::load(filename);
    } else if (suffixMatch(filename, ".pba")) {
        return FilePBA::load(filename);
    }

    panic("Unknown file format %s\n", filename.c_str());

    // keep compiler happy
    return Image();
}

void LoadFrames::help() {
    printf("\n-loadframes accepts a sequence of images and loads them as the frames of a\n"
           "single stack entry. See the help for -load for details on file formats.\n\n"
           "-loadframes cannot be used on raw float files. To achieve the same effect, cat\n"
           "the files together and load them as a single multi-frame image.\n\n"
           "Usage: ImageStack -loadframes foo*.jpg bar*.png\n\n");
}

bool LoadFrames::test() {
    Image a(123, 234, 3, 3);
    Noise::apply(a, 0, 1);
    FastBlur::apply(a, 1, 1, 1);

    string prefix = "_test";
    TempFile f1(prefix + "0.jpg");
    TempFile f2(prefix + "1.jpg");
    TempFile f3(prefix + "2.jpg");
    SaveFrames::apply(a, prefix + "%d.jpg", "100");
    vector<string> filenames;
    filenames.push_back(f1.name);
    filenames.push_back(f2.name);
    filenames.push_back(f3.name);
    Image b = LoadFrames::apply(filenames);
    return nearlyEqual(a, b);
}

void LoadFrames::parse(vector<string> args) {
    push(apply(args));
}


Image LoadFrames::apply(vector<string> args) {
    assert(args.size() > 0, "-loadframes requires at least one file argument.\n");

    Image im = Load::apply(args[0]);
    assert(im.frames == 1, "-loadframes can only load many single frame images\n");
    Image result(im.width, im.height, (int)args.size(), im.channels);
    result.frame(0).set(im);

    for (size_t i = 1; i < args.size(); i++) {
        im = Load::apply(args[i]);
        // check dimensions and channels match
        assert(im.frames == 1, "-loadframes can only load many single frame images\n");
        assert((im.width == result.width) &&
               (im.height == result.height) &&
               (im.channels == result.channels),
               "-loadframes can only load file sequences of matching width, height, and channel count\n");
        result.frame(i).set(im);
    }

    return result;
}


void LoadChannels::help() {
    pprintf("-loadchannels accepts a sequence of images and loads them as the channels of a"
            "single stack entry. See the help for -load for details on file formats.\n\n"
            "Usage: ImageStack -loadchannels foo*.jpg bar*.png\n\n");
}

bool LoadChannels::test() {
    Image a(123, 234, 1, 3);
    Noise::apply(a, 0, 1);
    FastBlur::apply(a, 1, 1, 1);

    string prefix = "_test";
    TempFile t1(prefix + "0.jpg");
    TempFile t2(prefix + "1.jpg");
    TempFile t3(prefix + "2.jpg");
    SaveChannels::apply(a, prefix + "%d.jpg", "100");
    vector<string> filenames;
    filenames.push_back(t1.name);
    filenames.push_back(t2.name);
    filenames.push_back(t3.name);
    Image b = LoadChannels::apply(filenames);
    return nearlyEqual(a, b);
}

void LoadChannels::parse(vector<string> args) {
    push(apply(args));
}


Image LoadChannels::apply(vector<string> args) {
    assert(args.size() > 0, "-loadchannels requires at least one file argument.\n");

    Image im = Load::apply(args[0]);
    assert(im.channels == 1, "-loadchannels can only load many single channel images\n");
    Image result(im.width, im.height, im.frames, (int)args.size());
    result.channel(0).set(im);

    for (size_t i = 1; i < args.size(); i++) {
        im = Load::apply(args[i]);
        // check dimensions and channels match
        assert(im.channels == 1, "-loadchannels can only load many single frame images\n");
        assert((im.width == result.width) &&
               (im.height == result.height) &&
               (im.frames == result.frames),
               "-loadchannels can only load file sequences of matching size\n");
        result.channel(i).set(im);
    }

    return result;
}


void Save::help() {
    printf("\n-save stores the image at the top of the stack to a file. The stack is not\n"
           "altered. The following file formats are supported:\n\n");

    printf("\n");
    FileTMP::help();
    printf("\n");
    FileHDR::help();
    printf("\n");
    FileJPG::help();
    printf("\n");
    FileEXR::help();
    printf("\n");
    FilePNG::help();
    printf("\n");
    FilePPM::help();
    printf("\n");
    FileTGA::help();
    printf("\n");
    FileWAV::help();
    printf("\n");
    FileTIFF::help();
    printf("\n");
    FileFLO::help();
    printf("\n");
    FileCSV::help();
    printf("\n");
    FilePBA::help();
    printf("\n");

    printf("Usage: ImageStack -load in.ppm -save out.jpg 98\n"
           "       ImageStack -load in.ppm -save out.jpg\n"
           "       ImageStack -load in.ppm -save out.ppm 16\n\n");

}

bool Save::test() {
    // tested by load
    return true;
}

void Save::parse(vector<string> args) {
    assert(args.size() == 1 || args.size() == 2, "-save requires exactly one or two arguments\n");
    if (args.size() == 1) { apply(stack(0), args[0], ""); }
    else if (args.size() == 2) { apply(stack(0), args[0], args[1]); }
}


void Save::apply(Image im, string filename, string arg) {
    if (suffixMatch(filename, ".tmp")) {
        if (arg == "") { arg = "float32"; }
        FileTMP::save(im, filename, arg);
    } else if (suffixMatch(filename, ".hdr")) {
        FileHDR::save(im, filename);
    } else if (suffixMatch(filename, ".jpg") ||
               suffixMatch(filename, ".jpeg")) {
        int quality;
        if (arg == "") quality = 90;
        else quality = readInt(arg);
        FileJPG::save(im, filename, quality);
    } else if (suffixMatch(filename, ".exr")) {
        string compression;
        if (arg == "") { compression = "piz"; }
        else { compression = arg; }
        FileEXR::save(im, filename, compression);
    } else if (suffixMatch(filename, ".png")) {
        FilePNG::save(im, filename);
    } else if (suffixMatch(filename, ".tga")) {
        FileTGA::save(im, filename);
    } else if (suffixMatch(filename, ".wav")) {
        FileWAV::save(im, filename);
    } else if (suffixMatch(filename, ".ppm") ||
               suffixMatch(filename, ".pgm")) {
        int bitdepth;
        if (arg == "") bitdepth = 16;
        else bitdepth = readInt(arg);
        FilePPM::save(im, filename, bitdepth);
    } else if (suffixMatch(filename, ".tiff") ||
               suffixMatch(filename, ".tif")) {
        FileTIFF::save(im, filename, arg);
    } else if (suffixMatch(filename, ".flo")) {
        FileFLO::save(im, filename);
    } else if (suffixMatch(filename, ".csv")) {
        FileCSV::save(im, filename);
    } else if (suffixMatch(filename, ".pba")) {
        FilePBA::save(im, filename);
    } else {
        panic("Unknown file format %s\n", filename.c_str());
    }
}

void SaveFrames::help() {
    printf("\n-saveframes takes a printf style format argument, and saves all the frames in\n"
           "the current image as separate files. See the help for save for details on file\n"
           "formats.\n\n"
           "Usage: ImageStack -loadframes *.jpg -saveframes frame%%03d.png\n"
           "       ImageStack -loadframes *.jpg -saveframes frame%%03d.ppm 16\n\n");
}

bool SaveFrames::test() {
    // Tested by LoadFrames
    return true;
}

void SaveFrames::parse(vector<string> args) {
    assert(args.size() == 1 || args.size() == 2, "-saveframes takes one or two arguments.\n");
    if (args.size() == 1) {
        apply(stack(0), args[0], "");
    } else {
        apply(stack(0), args[0], args[1]);
    }
}

// Microsoft has an underscore in front of snprintf for some reason
#ifdef _MSC_VER
#ifndef snprintf
#define snprintf _snprintf
#endif
#endif

void SaveFrames::apply(Image im, string pattern, string arg) {
    char filename[4096];

    for (int t = 0; t < im.frames; t++) {
        snprintf(filename, 4096, pattern.c_str(), t);
        Save::apply(im.frame(t), filename, arg);
    }
}

void SaveChannels::help() {
    pprintf("-savechannels takes a printf style format argument, and saves all the"
            " channels in the current image as separate files. See the help for"
            " save for details on file formats.\n"
            "\n"
            "Usage: ImageStack -loadchannels *.jpg -savechannels frame%%03d.png\n"
            "       ImageStack -loadchannels *.jpg -savechannels frame%%03d.ppm 16\n\n");
}

bool SaveChannels::test() {
    // Just calls save repeatedly
    return true;
}

void SaveChannels::parse(vector<string> args) {
    assert(args.size() == 1 || args.size() == 2, "-savechannels takes one or two arguments.\n");
    if (args.size() == 1) {
        apply(stack(0), args[0], "");
    } else {
        apply(stack(0), args[0], args[1]);
    }
}

void SaveChannels::apply(Image im, string pattern, string arg) {
    char filename[4096];

    for (int c = 0; c < im.channels; c++) {
        snprintf(filename, 4096, pattern.c_str(), c);
        Save::apply(im.channel(c), filename, arg);
    }
}

void LoadBlock::help() {
    pprintf("-loadblock loads a rectangular portion of a .tmp file. It is roughly"
            " equivalent to a load followed by a crop, except that the file need not"
            " fit in memory. The nine arguments are filename, x, y, t, and c offsets"
            " within the file, then width, height, frames, and channels. If seven"
            " arguments are given, all channels are loaded. If five arguments are"
            " given, all frames are used and the arguments specify x and y. If three"
            " arguments are given, they specify frames and all x, y, and channels"
            " are loaded. Loading out of bounds from the tmp file is"
            " permitted. Undefined areas will be zero-filled.\n\n"
            "This example multiplies a 512x512x128x3 volume by two, without ever loading it\n"
            "all into memory:\n"
            "ImageStack -loadblock foo.tmp 0 0 0 0 512 512 64 3 \\\n"
            "           -scale 2 -saveblock foo.tmp 0 0 0 0\n"
            "ImageStack -loadblock foo.tmp 0 0 0 64 512 512 64 3 \\\n"
            "           -scale 2 -saveblock foo.tmp 0 0 64 0\n\n"
           );
}

bool LoadBlock::test() {

    Image a(123, 234, 3, 3);
    Noise::apply(a, 0, 1);
    FastBlur::apply(a, 1, 1, 1);

    TempFile f(string("_test") + ".tmp");

    CreateTmp::apply(f.name, 234, 342, 5, 5);

    SaveBlock::apply(a, f.name, 4, 3, 2, 1);

    // Check the saved region is correct
    Image b = LoadBlock::apply(f.name, 4, 3, 2, 1, 123, 234, 3, 3);
    if (!nearlyEqual(a, b)) return false;

    // Check other regions are zero
    b = LoadBlock::apply(f.name, 130, 0, 0, 0, 50, 50, 5, 5);
    Stats s(b);
    return (s.mean() == 0 && s.variance() == 0);
}

void LoadBlock::parse(vector<string> args) {
    int t = 0, x = 0, y = 0, c = 0;
    int frames = -1, width = -1, height = -1, channels = -1;
    if (args.size() == 3) {
        t = readInt(args[1]);
        frames = readInt(args[2]);
    } else if (args.size() == 5) {
        x = readInt(args[1]);
        y = readInt(args[2]);
        width = readInt(args[3]);
        height = readInt(args[4]);
    } else if (args.size() == 7) {
        x = readInt(args[1]);
        y = readInt(args[2]);
        t = readInt(args[3]);
        width = readInt(args[4]);
        height = readInt(args[5]);
        frames = readInt(args[6]);
    } else if (args.size() == 9) {
        x = readInt(args[1]);
        y = readInt(args[2]);
        t = readInt(args[3]);
        c = readInt(args[4]);
        width = readInt(args[5]);
        height = readInt(args[6]);
        frames = readInt(args[7]);
        channels = readInt(args[8]);
    } else {
        panic("-loadblock takes 3, 5, 7, or 9 arguments\n");
    }
    push(LoadBlock::apply(args[0], x, y, t, c, width, height, frames, channels));
}

// 64 bit systems may not have fseeko. In this case fseek is just fine.
#ifndef fseeko
#define fseeko fseek
#endif

// A wrapper around fread that raises an exception if it runs out of data
void fread_(void *ptr, size_t size, size_t n, FILE *f) {
    assert(fread(ptr, size, n, f) == n,
           "Unexpected end of file\n");
}

Image LoadBlock::apply(string filename, int xoff, int yoff, int toff, int coff,
                       int width, int height, int frames, int channels) {
    // peek in the header

    struct {
        int width, height, frames, channels, type;
    } header;
    FILE *f = fopen(filename.c_str(), "rb");
    assert(f, "Could not open file: %s", filename.c_str());
    fread_(&header, sizeof(int), 5, f);

    if (width    <= 0) { width    = header.width; }
    if (height   <= 0) { height   = header.height; }
    if (frames   <= 0) { frames   = header.frames; }
    if (channels <= 0) { channels = header.channels; }

    if (header.type != 0) {
        fclose(f);
        panic("-loadblock can only handle tmp files containing floating point data.\n");
        return Image();
    }

    // sanity check the header
    if (header.width < 1 ||
        header.height < 1 ||
        header.frames < 1 ||
        header.channels < 1) {
        fclose(f);
        panic("According the header of the tmp file, the image has dimensions %ix%ix%ix%i. "
              "Perhaps this is not a tmp file?\n", header.width, header.height, header.frames, header.channels);
        return Image();
    }

    Image out(width, height, frames, channels);

    off_t xStride = sizeof(float);
    off_t yStride = header.width*xStride;
    off_t tStride = header.height*yStride;
    off_t cStride = header.frames*tStride;
    const off_t headerBytes = 5*sizeof(int);

    int xmin = max(xoff, 0);
    int xmax = min(xoff+width, header.width);
    int cmin = max(coff, 0);
    int cmax = min(coff+channels, header.channels);
    int ymin = max(yoff, 0);
    int ymax = min(yoff+height, header.height);
    int tmin = max(toff, 0);
    int tmax = min(toff+frames, header.frames);

    for (int c = cmin; c < cmax; c++) {
        for (int t = tmin; t < tmax; t++) {
            for (int y = ymin; y < ymax; y++) {
                off_t offset = c*cStride + t*tStride + y*yStride + xmin*xStride + headerBytes;
                fseeko(f, offset, SEEK_SET);
                fread_(&out(xmin-xoff, y-yoff, t-toff, c-coff), sizeof(float), (xmax-xmin), f);
            }
        }
    }

    fclose(f);

    return out;
}


void SaveBlock::help() {
    pprintf("-saveblock overwrites a rectangular subblock of a .tmp file with the"
            " top of the stack. It is logically similar to a load, paste, save"
            " combination, but never loads the full tmp file. The five arguments are"
            " the filename, followed by the offset at which to paste the volume in"
            " x, y, t, and c. When given four arguments c is assumed to be"
            " zero. When given three arguments t is also set to zero. With two"
            " arguments, x, y, and c are set to zero.\n\n"
            "This example multiplies a 128x512x512x3 volume by two, without ever loading it\n"
            "all into memory:\n"
            "ImageStack -loadblock foo.tmp 0 0 0 0 64 512 512 3 \\\n"
            "           -scale 2 -saveblock foo.tmp 0 0 0 0\n"
            "ImageStack -loadblock foo.tmp 64 0 0 0 64 512 512 3 \\\n"
            "           -scale 2 -saveblock foo.tmp 0 0 0 0\n\n");
}

bool SaveBlock::test() {
    // tested by loadblock
    return true;
}

void SaveBlock::parse(vector<string> args) {
    int t = 0, x = 0, y = 0, c = 0;
    if (args.size() == 2) {
        t = readInt(args[1]);
    } else if (args.size() == 3) {
        x = readInt(args[1]);
        y = readInt(args[2]);
    } else if (args.size() == 4) {
        x = readInt(args[1]);
        y = readInt(args[2]);
        t = readInt(args[3]);
    } else if (args.size() == 5) {
        x = readInt(args[1]);
        y = readInt(args[2]);
        t = readInt(args[3]);
        c = readInt(args[4]);
    } else {
        panic("-saveblock takes 2 to 5 arguments\n");
    }

    SaveBlock::apply(stack(0), args[0], x, y, t, c);
}

void SaveBlock::apply(Image im, string filename, int xoff, int yoff, int toff, int coff) {
    // Peek in the header
    struct {
        int width, height, frames, channels, type;
    } header;
    FILE *f = fopen(filename.c_str(), "rb+");
    assert(f, "Could not open file: %s", filename.c_str());
    fread_(&header, sizeof(int), 5, f);

    if (header.type != 0) {
        fclose(f);
        panic("-loadblock can only handle tmp files containing floating point data.\n");
        return;
    }

    // sanity check the header
    if (header.width < 1 || header.height < 1 || header.frames < 1 || header.channels < 1) {
        fclose(f);
        panic("According the header of the tmp file, the image has dimensions %ix%ix%ix%i. "
              "Perhaps this is not a tmp file?\n", header.width, header.height, header.frames, header.channels);
        return;
    }

    off_t xStride = sizeof(float);
    off_t yStride = header.width*xStride;
    off_t tStride = header.height*yStride;
    off_t cStride = header.frames*tStride;
    const int headerBytes = 5*sizeof(float);

    int xmin = max(xoff, 0);
    int xmax = min(xoff+im.width, header.width);
    int cmin = max(coff, 0);
    int cmax = min(coff+im.channels, header.channels);
    int ymin = max(yoff, 0);
    int ymax = min(yoff+im.height, header.height);
    int tmin = max(toff, 0);
    int tmax = min(toff+im.frames, header.frames);

    for (int c = cmin; c < cmax; c++) {
        for (int t = tmin; t < tmax; t++) {
            for (int y = ymin; y < ymax; y++) {
                off_t offset = c*cStride + t*tStride + y*yStride + xmin*xStride + headerBytes;
                fseeko(f, offset, SEEK_SET);
                fwrite(&im(xmin-xoff, y-yoff, t-toff, c-coff), sizeof(float), (xmax-xmin), f);
            }
        }

    }

    fclose(f);
}


void CreateTmp::help() {
    pprintf("-createtmp creates a zero filled floating point .tmp file of the specified"
            " dimensions. It can be used to create tmp files larger than can fit in"
            " memory. The five arguments are the filename, width, height, frames and"
            " channels. If only four arguments are specified, frames is assumed to"
            " be one.\n\n"
            "The following example creates a giant volume, and fills some of it with noise:\n"
            "ImageStack -createtmp volume.tmp 1024 1024 1024 1 \\\n"
            "           -push 256 256 256 1 -noise \\\n"
            "           -saveblock volume.tmp 512 512 512 0 \n\n");
}

bool CreateTmp::test() {
    // Tested by LoadBlock
    return true;
}

void CreateTmp::parse(vector<string> args) {
    int frames = 1, width = 1, height = 1, channels = 1;
    if (args.size() == 5) {
        width = readInt(args[1]);
        height = readInt(args[2]);
        frames = readInt(args[3]);
        channels = readInt(args[4]);
    } else if (args.size() == 4) {
        width = readInt(args[1]);
        height = readInt(args[2]);
        channels = readInt(args[3]);
    } else {
        panic("-createtmp takes either four or five arguments\n");
    }

    CreateTmp::apply(args[0], width, height, frames, channels);
}


void CreateTmp::apply(string filename, int width, int height, int frames, int channels) {

    assert(frames > 0 && width > 0 && height > 0 && channels > 0,
           "Some of the specified dimensions are less than 1\n");


    FILE *f = fopen(filename.c_str(), "wb");
    assert(f, "Could not open/create file %s\n", filename.c_str());

    int header[] = {width, height, frames, channels, 0};
    fwrite(header, sizeof(float), 5, f);
    vector<float> scanline(width, 0);

    for (int i = 0; i < height*frames*channels; i++) {
        fwrite(&scanline[0], sizeof(float), width, f);
    }

    fclose(f);
}

void LoadArray::help() {
    pprintf("-loadarray loads raw arrays of various data types. It takes 6"
            " arguments. The first is the filename to be loaded, the second is the"
            " data type, which must be one of: int8, uint8, int16, uint16, int32,"
            " uint32, float32, float64, or equivalently: char, unsigned char, short,"
            " unsigned short, int, unsigned int, float, double. The last four"
            " arguments specify the dimensions in the order width, height, frames,"
            " and channels. Bear in mind that ImageStack stores values internally as"
            " 32 bits floats, so information will be lost when double arrays are"
            " loaded. Integer formats are not scaled to lie between zero and one,"
            " this must be done manually with the -scale operation.\n"
            "\n"
            "Usage: ImageStack -loadarray foo.bar uint8 640 480 1 3\n\n");
}

namespace {
template<typename T>
bool testLoadArray() {
    Image a(123, 234, 3, 4);
    Noise::apply(a, 0, 12324);
    TempFile f;
    SaveArray::apply<T>(a, f.name);
    Image b = LoadArray::apply<T>(f.name, 123, 234, 3, 4);

    for (int i = 0; i < 100; i++) {
        int x = randomInt(0, a.width-1);
        int y = randomInt(0, a.height-1);
        int t = randomInt(0, a.frames-1);
        int c = randomInt(0, a.channels-1);
        // LoadArray/SaveArray should be exactly equivalent to a C-style cast
        float saved = (float)((T)(a(x, y, t, c)));
        float loaded = b(x, y, t, c);
        if (loaded != saved) {
            printf("%f vs %f\n", saved, loaded);
            return false;
        }
    }
    return true;
}
};

bool LoadArray::test() {

    if (!testLoadArray<uint32_t>()) return false;
    if (!testLoadArray<int32_t>()) return false;
    if (!testLoadArray<uint16_t>()) return false;
    if (!testLoadArray<int16_t>()) return false;
    if (!testLoadArray<uint8_t>()) return false;
    if (!testLoadArray<int8_t>()) return false;
    if (!testLoadArray<float>()) return false;
    if (!testLoadArray<double>()) return false;

    return true;
}

void LoadArray::parse(vector<string> args) {
    assert(args.size() == 6, "-loadarray takes 6 arguments\n");
    string filename = args[0];
    string type = args[1];
    int frames = readInt(args[4]), width = readInt(args[2]);
    int height = readInt(args[3]), channels = readInt(args[5]);
    if (type == "int8" || type == "char") {
        push(apply<int8_t>(filename, width, height, frames, channels));
    } else if (type == "uint8" || type == "unsigned char") {
        push(apply<uint8_t>(filename, width, height, frames, channels));
    } else if (type == "int16" || type == "short") {
        push(apply<int16_t>(filename, width, height, frames, channels));
    } else if (type == "uint16" || type == "unsigned short") {
        push(apply<uint16_t>(filename, width, height, frames, channels));
    } else if (type == "int32" || type == "int") {
        push(apply<int32_t>(filename, width, height, frames, channels));
    } else if (type == "uint32" || type == "unsigned int") {
        push(apply<uint32_t>(filename, width, height, frames, channels));
    } else if (type == "float32" || type == "float") {
        push(apply<float>(filename, width, height, frames, channels));
    } else if (type == "float64" || type == "double") {
        push(apply<double>(filename, width, height, frames, channels));
    } else {
        panic("Unknown type %s\n", type.c_str());
    }
}

template<typename T>
Image LoadArray::apply(string filename, int width, int height, int frames, int channels) {
    FILE *f = fopen(filename.c_str(), "rb");

    assert(f, "Could not open file %s", filename.c_str());

    Image im(width, height, frames, channels);

    vector<T> rawData(width);

    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                fread_(&rawData[0], sizeof(T), width, f);
                for (int x = 0; x < im.width; x++) {
                    im(x, y, t, c) = (float)rawData[x];
                }
            }
        }
    }

    fclose(f);

    return im;
}

void SaveArray::help() {
    printf("-savearray saves raw arrays of various data types. It takes 2 arguments. The first\n"
           "is the filename to be loaded, the second is the data type, which must be one of:\n"
           "int8, uint8, int16, uint16, int32, uint32, float32, float64, or equivalently:\n"
           "char, unsigned char, short, unsigned short, int, unsigned int, float, double.\n\n"
           "Bear in mind that ImageStack stores values internally as 32 bit floats, so\n"
           "saving in double format does not give you higher fidelity.\n\n"
           "Integer type formats are not scaled from the range zero to one. This must be\n"
           "done manually using the -scale operation.\n\n"
           "Usage: ImageStack -load in.jpg -savearray out.float float32\n\n");
}

bool SaveArray::test() {
    // Tested by LoadArray
    return true;
}

void SaveArray::parse(vector<string> args) {
    assert(args.size() == 2, "-savearray takes 2 arguments\n");
    string filename = args[0];
    string type = args[1];
    if (type == "int8" || type == "char") {
        apply<int8_t>(stack(0), filename);
    } else if (type == "uint8" || type == "unsigned char") {
        apply<uint8_t>(stack(0), filename);
    } else if (type == "int16" || type == "short") {
        apply<int16_t>(stack(0), filename);
    } else if (type == "uint16" || type == "unsigned short") {
        apply<uint16_t>(stack(0), filename);
    } else if (type == "int32" || type == "int") {
        apply<int32_t>(stack(0), filename);
    } else if (type == "uint32" || type == "unsigned int") {
        apply<uint32_t>(stack(0), filename);
    } else if (type == "float32" || type == "float") {
        apply<float>(stack(0), filename);
    } else if (type == "float64" || type == "double") {
        apply<double>(stack(0), filename);
    } else {
        panic("Unknown type %s\n", type.c_str());
    }
}


template<typename T>
void SaveArray::apply(Image im, string filename) {
    FILE *f = fopen(filename.c_str(), "wb");

    assert(f, "Could not open file %s\n", filename.c_str());

    vector<T> rawData(im.width);

    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width; x++) {
                    rawData[x] = (T)(im(x, y, t, c));
                }
                fwrite(&rawData[0], sizeof(T), im.width, f);
            }
        }
    }

    fclose(f);
}
#include "footer.h"
