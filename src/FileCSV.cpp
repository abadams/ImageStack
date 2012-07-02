#include "main.h"
#include "File.h"
#include "header.h"

namespace FileCSV {
void help() {
    pprintf(".csv files. These contain comma-separated floating point values in"
            " text. Each scanline of the image corresponds to a line in the file. x"
            " and c are thus conflated, as are y and t. When loading csv files,"
            " ImageStack assumes 1 channel and 1 frame.\n");
}


Image load(string filename) {
    // calculate the number of rows and columns in the file
    FILE *f = fopen(filename.c_str(), "r");

    // how many commas in the first line?
    int width = 1;
    int c, last;
    do {
        c = fgetc(f);
        if (c == ',') { width++; }
    } while (c != '\n' && c != EOF);

    // how many lines in the file?
    int height = 1;
    do {
        last = c;
        c = fgetc(f);
        if (c == '\n' && last != '\n') { height++; }
    } while (c != EOF);

    // go back to the start and start reading data
    fseek(f, 0, SEEK_SET);

    Image out(width, height, 1, 1);

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width-1; x++) {
            assert(fscanf(f, "%20f,", &out(x, y)) == 1, "Failed to parse file\n");
        }
        assert(fscanf(f, "%20f", &out(width-1, y)) == 1, "Failed to parse file\n");
    }

    fclose(f);

    return out;
}

void save(Image im, string filename) {
    FILE *f = fopen(filename.c_str(), "w");

    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width-1; x++) {
                    fprintf(f, "%10.10f, ", im(x, y, t, c));
                }
                fprintf(f, "%10.10f\n", im(im.width-1, y, t, c));
            }
        }
    }

    fclose(f);
}
}
#include "footer.h"
