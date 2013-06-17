#include "main.h"
#include "File.h"
namespace ImageStack {

#ifdef NO_PNG
namespace FilePNG {
#include "FileNotImplemented.h"
}
#else

namespace FilePNG {

#define PNG_DEBUG 3
#include <png.h>

void help() {
    pprintf(".png files. These have a bit depth of 8 or 16, and may have 1-4 channels. They may only have 1 frame.");
}

Image load(string filename) {
    png_byte header[8];        // 8 is the maximum size that can be checked

    /* open file and test for it being a png */
    FILE *f = fopen(filename.c_str(), "rb");
    assert(f, "File %s could not be opened for reading\n", filename.c_str());
    assert(fread(header, 1, 8, f) == 8, "File ended before end of header\n");
    assert(!png_sig_cmp(header, 0, 8), "File %s is not recognized as a PNG file\n", filename.c_str());

    /* initialize stuff */
    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    assert(png_ptr, "[read_png_file] png_create_read_struct failed\n");

    png_infop info_ptr = png_create_info_struct(png_ptr);
    assert(info_ptr, "[read_png_file] png_create_info_struct failed\n");

    assert(!setjmp(png_jmpbuf(png_ptr)), "[read_png_file] Error during init_io\n");

    png_init_io(png_ptr, f);
    png_set_sig_bytes(png_ptr, 8);

    png_read_info(png_ptr, info_ptr);

    int width = png_get_image_width(png_ptr, info_ptr);
    int height = png_get_image_height(png_ptr, info_ptr);
    int channels = png_get_channels(png_ptr, info_ptr);
    int bit_depth = png_get_bit_depth(png_ptr, info_ptr);

    // Expand low-bpp images to have only 1 pixel per byte (As opposed to tight packing)
    if (bit_depth < 8) {
        png_set_packing(png_ptr);
    }

    Image im(width, height, 1, channels);

    //number_of_passes = png_set_interlace_handling(png_ptr);
    png_read_update_info(png_ptr, info_ptr);

    // read the file
    assert(!setjmp(png_jmpbuf(png_ptr)), "[read_png_file] Error during read_image\n");

    std::vector<png_bytep> row_pointers(im.height);
    int row_bytes = png_get_rowbytes(png_ptr, info_ptr);
    std::vector<png_byte> data(row_bytes * im.height);
    for (int y = 0; y < im.height; y++) {
        row_pointers[y] = &data[y * row_bytes];
    }

    png_read_image(png_ptr, &row_pointers[0]);

    fclose(f);

    // convert the data to floats
    if (bit_depth <= 8) {
        int bit_scale = 8/bit_depth;
        for (int y = 0; y < im.height; y++) {
            png_bytep srcPtr = row_pointers[y];
            for (int x = 0; x < im.width; x++) {
                for (int c = 0; c < im.channels; c++) {
                    im(x, y, c) = LDRtoHDR(bit_scale* (*srcPtr++));
                }
            }
        }
    } else if (bit_depth == 16) {
        for (int y = 0; y < im.height; y++) {
            png_bytep srcPtr = row_pointers[y];
            for (int x = 0; x < im.width; x++) {
                for (int c = 0; c < im.channels; c++) {
                    unsigned short val = srcPtr[0]*256 + srcPtr[1];
                    im(x, y, c) = LDR16toHDR(val);
                    srcPtr += 2;
                }
            }
        }
    }

    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);

    return im;
}


  void save(Image im, string filename, int bits) {
    assert(bits == 8 || bits == 16, "Can only save 8 or 16 bit pngs\n");
    assert(im.frames == 1, "Can't save a multi-frame PNG image\n");
    assert(im.channels > 0 && im.channels < 5,
           "Imagestack can't write PNG files that have other than 1, 2, 3, or 4 channels\n");

    png_byte color_types[4] = {PNG_COLOR_TYPE_GRAY, PNG_COLOR_TYPE_GRAY_ALPHA,
                               PNG_COLOR_TYPE_RGB,  PNG_COLOR_TYPE_RGB_ALPHA
                              };
    png_byte color_type = color_types[im.channels - 1];

    // open file
    FILE *f = fopen(filename.c_str(), "wb");
    assert(f, "[write_png_file] File %s could not be opened for writing\n", filename.c_str());

    // initialize stuff
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    assert(png_ptr, "[write_png_file] png_create_write_struct failed\n");

    png_infop info_ptr = png_create_info_struct(png_ptr);
    assert(info_ptr, "[write_png_file] png_create_info_struct failed\n");

    assert(!setjmp(png_jmpbuf(png_ptr)), "[write_png_file] Error during init_io\n");

    png_init_io(png_ptr, f);

    // write header
    assert(!setjmp(png_jmpbuf(png_ptr)), "[write_png_file] Error during writing header\n");

    png_set_IHDR(png_ptr, info_ptr, im.width, im.height,
                 bits, color_type, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

    png_write_info(png_ptr, info_ptr);

    // convert the floats to bytes
    int row_bytes = png_get_rowbytes(png_ptr, info_ptr);
    std::vector<png_bytep> row_pointers(im.height);
    std::vector<png_byte> data(row_bytes * im.height);
    for (int y = 0; y < im.height; y++) {
        png_bytep dstPtr = row_pointers[y] = &data[y*row_bytes];
        if (bits == 8) {
            for (int x = 0; x < im.width; x++) {
                for (int c = 0; c < im.channels; c++) {                
                    *dstPtr++ = (png_byte)(HDRtoLDR(im(x, y, c)));
                }
            }
        } else if (bits == 16) {
            for (int x = 0; x < im.width; x++) {
                for (int c = 0; c < im.channels; c++) {                
                    unsigned short val = HDRtoLDR16(im(x, y, c));
                    *dstPtr++ = (png_byte)(val >> 8);
                    *dstPtr++ = (png_byte)(val & 0xff);
                }
            }            
        }
    }

    // write data
    assert(!setjmp(png_jmpbuf(png_ptr)), "[write_png_file] Error during writing bytes");

    png_write_image(png_ptr, &row_pointers[0]);

    // finish write
    assert(!setjmp(png_jmpbuf(png_ptr)), "[write_png_file] Error during end of write");

    png_write_end(png_ptr, NULL);

    fclose(f);

    png_destroy_write_struct(&png_ptr, &info_ptr);
}


}

#endif
}
