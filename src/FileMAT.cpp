#include "main.h"
#include "File.h"

#ifndef _WIN32
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

namespace ImageStack {

namespace FileMAT {

void help() {
    pprintf(".mat files in Matlab's level 5 format (as used by Matlab v6). "
            "When saving, an optional second argument specifies the format. This"
            " may be any of int8, uint8, int16, uint16, int32, uint32, int64,"
            " uint64, float32, float64, or correspondingly char, unsigned"
            " char, short, unsigned short, int, unsigned int, float, or"
            " double. The default is float32.\n");
}


template<typename T>
void saveData(FILE *f, Image im) {

    vector<T> buf(im.width);

    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                float *srcPtr = &im(0, y, t, c);
                for (int x = 0; x < im.width; x++) {
                    buf[x] = (T)srcPtr[x];
                }
                fwrite(&buf[0], sizeof(T), im.width, f);
            }
        }
    }
}

template<>
void saveData<float>(FILE *f, Image im) {

    vector<float> buf(im.width*im.channels);

    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                fwrite(&im(0, y, t, c), sizeof(float), im.width, f);
            }
        }
    }
}

enum MatlabTypeCode {
    miINT8 = 1,
    miUINT8 = 2,
    miINT16 = 3,
    miUINT16 = 4,
    miINT32 = 5,
    miUINT32 = 6,
    miSINGLE = 7,
    miDOUBLE = 9,
    miINT64 = 12,
    miUINT64 = 13,
    miMATRIX = 14,
    miCOMPRESSED = 15,
    miUTF8 = 16,
    miUTF16 = 17,
    miUTF32 = 18
};

enum MatlabClassCode {
    mxCHAR_CLASS = 3,
    mxDOUBLE_CLASS = 6,
    mxSINGLE_CLASS = 7,
    mxINT8_CLASS = 8,
    mxUINT8_CLASS = 9,
    mxINT16_CLASS = 10,
    mxUINT16_CLASS = 11,
    mxINT32_CLASS = 12,
    mxUINT32_CLASS = 13,
    mxINT64_CLASS = 14,
    mxUINT64_CLASS = 15
};

void save(Image im, string filename, string type) {
    FILE *f = fopen(filename.c_str(), "wb");
    assert(f, "Could not write output file %s\n", filename.c_str());

    int elem_size;
    uint8_t class_code;
    uint32_t type_code;

    if (type == "float" || type == "float32") {
        type_code = miSINGLE;
        class_code = mxSINGLE_CLASS;
        elem_size = 4;
    } else if (type == "double" || type == "float64") {
        type_code = miDOUBLE;
        class_code = mxDOUBLE_CLASS;
        elem_size = 4;
    } else if (type == "uint8" || type == "unsigned char") {
        type_code = miUINT8;
        class_code = mxUINT8_CLASS;
        elem_size = 1;
    } else if (type == "int8" || type == "char") {
        type_code = miINT8;
        class_code = mxINT8_CLASS;
        elem_size = 1;
    } else if (type == "uint16" || type == "unsigned short") {
        type_code = miUINT16;
        class_code = mxUINT16_CLASS;
        elem_size = 2;
    } else if (type == "int16" || type == "short") {
        type_code = miINT16;
        class_code = mxINT16_CLASS;
        elem_size = 2;
    } else if (type == "uint32" || type == "unsigned int") {
        type_code = miUINT32;
        class_code = mxUINT32_CLASS;
        elem_size = 4;
    } else if (type == "int32" || type == "int") {
        type_code = miINT32;
        class_code = mxINT32_CLASS;
        elem_size = 4;
    } else if (type == "uint64") {
        type_code = miUINT64;
        class_code = mxUINT64_CLASS;
        elem_size = 8;
    } else if (type == "int64") {
        type_code = miINT64;
        class_code = mxINT64_CLASS;
        elem_size = 8;
    }

    // Pick a name for the array
    size_t idx = filename.rfind(".");
    string name = filename.substr(0, idx);
    uint32_t name_size = (int)name.size();
    while (name.size() & 0x7) name += '\0';

    char header[128] = "MATLAB 5.0 MAT-file, produced by ImageStack";
    int len = strlen(header);
    memset(header + len, ' ', sizeof(header) - len);

    // Version
    *((uint16_t *)(header + 124)) = 0x0100;

    // Endianness check
    header[126] = 'I';
    header[127] = 'M';

    fwrite(header, sizeof(header), 1, f);

    uint32_t payload_bytes = im.width * im.height * im.frames * im.channels * elem_size;

    // Matrix header
    uint32_t matrix_header[2] = {
        miMATRIX, 16 + 24 + 8 + (uint32_t)name.size() + 8 + (uint32_t)payload_bytes
    };
    fwrite(matrix_header, sizeof(matrix_header), 1, f);

    // Array flags
    uint32_t flags[4] = {
        miUINT32, 8, class_code, 1
    };
    fwrite(flags, sizeof(flags), 1, f);

    // Shape
    int32_t shape[6] = {
        miINT32, 16, im.width, im.height, im.frames, im.channels
    };
    fwrite(&shape, sizeof(shape), 1, f);

    // Name
    uint32_t name_header[2] = {
        miINT8, name_size
    };

    fwrite(name_header, sizeof(name_header), 1, f);
    fwrite(&name[0], name.size(), 1, f);

    uint32_t padding_bytes = 7 - ((payload_bytes - 1) & 7);

    // Payload header
    uint32_t payload_header[2] = {
        type_code, payload_bytes
    };
    fwrite(payload_header, sizeof(payload_header), 1, f);

    // Payload
    switch (type_code) {
    case miSINGLE:
        saveData<float>(f, im);
        break;
    case miDOUBLE:
        saveData<double>(f, im);
        break;
    case miUINT8:
        saveData<uint8_t>(f, im);
        break;
    case miINT8:
        saveData<int8_t>(f, im);
        break;
    case miUINT16:
        saveData<uint16_t>(f, im);
        break;
    case miINT16:
        saveData<int16_t>(f, im);
        break;
    case miUINT32:
        saveData<uint32_t>(f, im);
        break;
    case miINT32:
        saveData<int32_t>(f, im);
        break;
    case miUINT64:
        saveData<uint64_t>(f, im);
        break;
    case miINT64:
        saveData<int64_t>(f, im);
        break;
    }

    // Padding
    assert(padding_bytes < 8, "Too much padding: %d\n", padding_bytes);
    uint64_t padding = 0;
    fwrite(&padding, padding_bytes, 1, f);

    fclose(f);
}

namespace {
void checked_fread(void *ptr, size_t size, size_t nmemb, FILE *stream) {
    size_t items_read = fread(ptr, size, nmemb, stream);
    size_t items_expected = nmemb;
    if (items_read != items_expected) {
        fprintf(stderr, "Unexpected end of file. Tried to read %llu, got %llu\n",
                (long long unsigned)(items_expected), (long long unsigned)items_read);
        abort();
    }
}
}

template<typename T>
Image loadData(FILE *f, int width, int height, int frames, int channels) {
    Image im(width, height, frames, channels);

    vector<T> scanline(width);
    for (int c = 0; c < channels; c++) {
        for (int t = 0; t < frames; t++) {
            for (int y = 0; y < height; y++) {
                checked_fread(&scanline[0], sizeof(T), width, f);
                for (int x = 0; x < width; x++) {
                    im(x, y, t, c) = (float)scanline[x];
                }
            }
        }
    }

    return im;
}

template<>
Image loadData<float>(FILE *f, int width, int height, int frames, int channels) {
    Image im(width, height, frames, channels);

    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames; t++) {
            for (int y = 0; y < im.height; y++) {
                checked_fread(&im(0, y, t, c), sizeof(float), im.width, f);
            }
        }
    }

    return im;
}

Image load(string filename) {
    static map<string, FILE *> open_fifos;

    FILE *f = NULL;

    map<string, FILE *>::iterator it = open_fifos.find(filename);
    if (it == open_fifos.end()) {
        f = fopen(filename.c_str(), "rb");
    } else {
        f = it->second;
    }

    assert(f, "Could not open file %s\n", filename.c_str());
    uint8_t header[128];
    checked_fread(&header, sizeof(header), 1, f);

    // Matrix header
    uint32_t matrix_header[2];
    checked_fread(matrix_header, sizeof(matrix_header), 1, f);
    assert(matrix_header[0] == miMATRIX, "%d) Failed to parse .mat file %s", __LINE__, filename.c_str());

    // Array flags
    uint32_t flags[4];
    checked_fread(flags, sizeof(flags), 1, f);
    assert(flags[0] == miUINT32, "%d) Failed to parse .mat file %s", __LINE__, filename.c_str());
    assert(flags[1] == 8, "%d) Failed to parse .mat file %s", __LINE__, filename.c_str());

    // Shape
    uint32_t shape_header[2];
    checked_fread(shape_header, sizeof(shape_header), 1, f);
    assert(shape_header[0] == miINT32, "%d) Failed to parse .mat file %s", __LINE__, filename.c_str());
    int dims = (shape_header[1])/4;
    assert(dims <= 4, "ImageStack cannot load .mat files with more than four dimensions");
    uint32_t shape[4] = {1, 1, 1, 1};
    checked_fread(&shape, 4, dims, f);

    // Skip over the name
    uint32_t name_header[2];
    checked_fread(&name_header, sizeof(name_header), 1, f);
    if (name_header[0] >> 16) {
        // Name must be fewer than 4 chars, and so the whole name
        // field was stored packed into 8 bytes
    } else {
        assert(name_header[0] == miINT8, "%d) Failed to parse .mat file %s", __LINE__, filename.c_str());
        int name_bytes = name_header[1];
        while (name_bytes & 7) name_bytes++;
        fseek(f, name_bytes, SEEK_CUR);
    }

    // Payload header
    uint32_t payload_header[2];
    checked_fread(payload_header, sizeof(payload_header), 1, f);
    uint32_t type_code = payload_header[0];
    Image im;

    switch (type_code) {
    case miSINGLE:
        im = loadData<float>(f, shape[0], shape[1], shape[2], shape[3]);
        break;
    case miDOUBLE:
        im = loadData<double>(f, shape[0], shape[1], shape[2], shape[3]);
        break;
    case miINT8:
        im = loadData<int8_t>(f, shape[0], shape[1], shape[2], shape[3]);
        break;
    case miUINT8:
        im = loadData<uint8_t>(f, shape[0], shape[1], shape[2], shape[3]);
        break;
    case miINT16:
        im = loadData<int16_t>(f, shape[0], shape[1], shape[2], shape[3]);
        break;
    case miUINT16:
        im = loadData<uint16_t>(f, shape[0], shape[1], shape[2], shape[3]);
        break;
    case miINT32:
        im = loadData<int32_t>(f, shape[0], shape[1], shape[2], shape[3]);
        break;
    case miUINT32:
        im = loadData<uint32_t>(f, shape[0], shape[1], shape[2], shape[3]);
        break;
    case miINT64:
        im = loadData<int64_t>(f, shape[0], shape[1], shape[2], shape[3]);
        break;
    case miUINT64:
        im = loadData<uint64_t>(f, shape[0], shape[1], shape[2], shape[3]);
        break;
    default:
        printf("Unknown type code %d.\n", type_code);
    }

    bool is_fifo = false;

#ifndef _WIN32
    {
        struct stat buf {0};
        is_fifo = (fstat(fileno(f), &buf) == 0) && (S_ISFIFO(buf.st_mode));
    }
#endif

    // If the file is actually a fifo, keep it open so that we don't
    // crash the writer.
    if (is_fifo) {
        open_fifos[filename] = f;
    } else {
        fclose(f);
    }

    return im;
}
}
}
