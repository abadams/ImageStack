#include "main.h"
#include "File.h"
#include "header.h"

namespace FilePBA {

void help() {
    pprintf(".pba files. This format is a human-readable space-separated 2D array of"
            " numbers. It is used by petabricks. Width becomes columns"
            " of the file, and height, frames, and channels become rows.\n");
}

void save(Image im, string filename) {
    FILE *f = fopen(filename.c_str(), "w");
    assert(f, "Could not write output file %s\n", filename.c_str());
    // write the dimensions
    fprintf(f, "SIZE");
    if (im.channels != 1) { fprintf(f, " %d", im.channels); }
    if (im.width != 1) { fprintf(f, " %d", im.width); }
    if (im.height != 1) { fprintf(f, " %d", im.height); }
    if (im.frames != 1) { fprintf(f, " %d", im.frames); }
    fprintf(f, "\n");

    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width; x++) {
                    fprintf(f, "%f ", im(x, y, t, c));
                }
                fprintf(f, "\n");
            }
            fprintf(f, "\n");
        }
        fprintf(f, "\n");
    }

    fclose(f);
}

Image load(string filename) {
    FILE *f = fopen(filename.c_str(), "r");
    assert(f, "Could not read input file %s\n", filename.c_str());

    // read the dimensions
    char header[256];
    assert(fgets(header, 256, f) == header,
           "Could not read header of %s\n", filename.c_str());

    int ch, wi, he, fr;
    int ret = sscanf(header, "SIZE %20d %20d %20d %20d", &ch, &wi, &he, &fr);
    if (ret == 0) { // 0D
        ch = wi = he = fr = 1;
    } else if (ret == 1) { // 1D
        wi = ch;
        ch = he = fr = 1;
    } else if (ret == 2) { // 2D
        he = wi;
        wi = ch;
        ch = fr = 1;
    } else if (ret == 3) {
        fr = 1;
    } else if (ret == 4) {

    } else {
        panic("Could not parse header of %s\n", filename.c_str());
    }



    Image im(wi, he, fr, ch);

    // read the data
    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width; x++) {
                    assert(fscanf(f, "%20f", &im(x, y, t, c)) == 1,
                           "Unexpected end of file reading %s\n", filename.c_str());
                }
            }
        }
    }

    fclose(f);

    return im;
}

}
#include "footer.h"
