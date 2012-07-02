#include "main.h"
#include "Geometry.h"
#include "Stack.h"
#include "Arithmetic.h"
#include "Statistics.h"
#include "header.h"

void Upsample::help() {
    pprintf("-upsample multiplies the width, height, and frames of the current"
            " image by the given integer arguments. It uses nearest neighbor"
            " interpolation. For a slower, high-quality resampling method, use"
            " -resample instead.\n\n"
            "-upsample x y is interpreted as -upsample x y 1\n"
            "-upsample x is interpreted as -upsample x x 1\n"
            "-upsample is interpreted as -upsample 2 2 1\n\n"
            "Usage: ImageStack -load a.tga -upsample 3 2 -save b.tga\n\n");
}

bool Upsample::test() {
    Image a(12, 23, 34, 2);
    Noise::apply(a, -2, 3);
    Image b = Upsample::apply(a, 2, 3, 4);
    if (b.width != a.width*2) return false;
    if (b.height != a.height*3) return false;
    if (b.frames != a.frames*4) return false;
    if (b.channels != a.channels) return false;
    for (int i = 0; i < 10; i++) {
        int x = randomInt(0, b.width-1);
        int y = randomInt(0, b.height-1);
        int t = randomInt(0, b.frames-1);
        int c = randomInt(0, a.channels-1);
        if (b(x, y, t, c) != a(x/2, y/3, t/4, c)) return false;
    }
    return true;
}

void Upsample::parse(vector<string> args) {
    int boxWidth = 2, boxHeight = 2, boxFrames = 1;
    assert(args.size() <= 3, "-upsample takes three or fewer arguments\n");
    if (args.size() == 3) {
        boxWidth = readInt(args[0]);
        boxHeight = readInt(args[1]);
        boxFrames = readInt(args[2]);
    } else if (args.size() == 2) {
        boxWidth = readInt(args[0]);
        boxHeight = readInt(args[1]);
    } else if (args.size() == 1) {
        boxWidth = boxHeight = readInt(args[0]);
    }

    Image im = apply(stack(0), boxWidth, boxHeight, boxFrames);
    pop();
    push(im);
}

Image Upsample::apply(Image im, int boxWidth, int boxHeight, int boxFrames) {

    Image out(im.width*boxWidth, im.height*boxHeight, im.frames*boxFrames, im.channels);

    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < out.frames; t++) {
            int it = t / boxFrames;
            for (int y = 0; y < out.height; y++) {
                int iy = y / boxHeight;
                for (int x = 0; x < out.width; x++) {
                    int ix = x / boxWidth;
                    out(x, y, t, c) = im(ix, iy, it, c);
                }
            }
        }
    }

    return out;

}

void Downsample::help() {
    pprintf("-downsample divides the width, height, and frames of the current image"
            " by the given integer arguments. It averages rectangles to get the new"
            " values.\n\n"
            "-downsample x y is interpreted as -downsample x y 1\n"
            "-downsample x is interpreted as -downsample x x 1\n"
            "-downsample is interpreted as -downsample 2 2 1\n\n"
            "Usage: ImageStack -load a.tga -downsample 3 2 -save b.tga\n\n");
}

bool Downsample::test() {
    Image a(12, 23, 34, 2);
    Noise::apply(a, -2, 3);
    Image b = Upsample::apply(a, 2, 3, 4);
    Image a2 = Downsample::apply(b, 2, 3, 4);
    if (a.width != a2.width) return false;
    if (a.height != a2.height) return false;
    if (a.frames != a2.frames) return false;
    if (a.channels != a2.channels) return false;
    for (int i = 0; i < 10; i++) {
        int x = randomInt(0, a2.width-1);
        int y = randomInt(0, a2.height-1);
        int t = randomInt(0, a2.frames-1);
        int c = randomInt(0, a2.channels-1);
        if (!nearlyEqual(a(x, y, t, c), a2(x, y, t, c))) return false;
    }
    return true;
}

void Downsample::parse(vector<string> args) {
    int boxWidth = 2, boxHeight = 2, boxFrames = 1;
    assert(args.size() <= 3, "-downsample takes three or fewer arguments\n");
    if (args.size() == 3) {
        boxWidth = readInt(args[0]);
        boxHeight = readInt(args[1]);
        boxFrames = readInt(args[2]);
    } else if (args.size() == 2) {
        boxWidth = readInt(args[0]);
        boxHeight = readInt(args[1]);
    } else if (args.size() == 1) {
        boxWidth = boxHeight = readInt(args[0]);
    }

    Image im = apply(stack(0), boxWidth, boxHeight, boxFrames);
    pop();
    push(im);
}

Image Downsample::apply(Image im, int boxWidth, int boxHeight, int boxFrames) {

    //if (!((im.width % boxWidth == 0) && (im.height % boxHeight == 0) && (im.frames % boxFrames == 0))) {
    //printf("Warning: Image dimensions are not a multiple of the downsample size. Ignoring some pixels.\n");
    //}

    int newWidth = im.width / boxWidth;
    int newHeight = im.height / boxHeight;
    int newFrames = im.frames / boxFrames;
    float scale = 1.0f / (boxWidth * boxHeight * boxFrames);

    Image out(newWidth, newHeight, newFrames, im.channels);

    for (int c = 0; c < out.channels; c++) {
        for (int t = 0; t < out.frames; t++) {
            for (int y = 0; y < out.height; y++) {
                for (int x = 0; x < out.width; x++) {
                    float val = 0;
                    for (int dt = 0; dt < boxFrames; dt++) {
                        for (int dy = 0; dy < boxHeight; dy++) {
                            for (int dx = 0; dx < boxWidth; dx++) {
                                val += im(x*boxWidth+dx, y*boxHeight+dy, t*boxFrames+dt, c);
                            }
                        }
                    }
                    out(x, y, t, c) = val * scale;
                }
            }
        }
    }

    return out;
}

void Resample::help() {
    printf("-resample resamples the input using a 3-lobed Lanczos filter. When"
           " given three arguments, it produces a new volume of the given width,"
           " height, and frames. When given two arguments, it produces a new volume"
           " of the given width and height, with the same number of frames.\n\n"
           "Usage: ImageStack -loadframes f*.tga -resample 20 50 50 -saveframes f%%03d.tga\n\n");
}

bool Resample::test() {
    Image a(4, 8, 3, 3);
    Noise::apply(a, 0, 1);
    a = Resample::apply(a, 50, 50, 3);
    Image b = Resample::apply(a, 225, 175, 3);
    if (b.width != 225 || b.height != 175 || b.frames != 3) return false;
    Image a2 = Resample::apply(b, 50, 50, 3);
    Stats s1(a), s2(a2);
    if (!nearlyEqual(s1.mean(), s2.mean())) return false;
    if (!nearlyEqual(s1.variance(), s2.variance())) return false;

    for (int c = 0; c < a.channels; c++) {
        for (int t = 0; t < a.frames; t++) {
            for (int y = 6; y < a.height-7; y++) {
                for (int x = 6; x < a.width-7; x++) {
                    if (fabs(a(x, y, t, c) - a2(x, y, t, c)) > 0.01) return false;
                }
            }
        }
    }
    return true;
}

void Resample::parse(vector<string> args) {

    if (args.size() == 2) {
        Image im = apply(stack(0), readInt(args[0]), readInt(args[1]));
        pop();
        push(im);
    } else if (args.size() == 3) {
        Image im = apply(stack(0), readInt(args[0]), readInt(args[1]), readInt(args[2]));
        pop();
        push(im);
    } else {
        panic("-resample takes two or three arguments\n");
    }

}

Image Resample::apply(Image im, int width, int height) {
    if (height != im.height && width != im.width) {
        Image tmp = resampleY(im, height);
        return resampleX(tmp, width);
    } else if (width != im.width) {
        return resampleX(im, width);
    } else if (height != im.height) {
        return resampleY(im, height);
    }
    return im;
}

Image Resample::apply(Image im, int width, int height, int frames) {
    if (frames != im.frames) {
        Image tmp = resampleT(im, frames);
        return apply(tmp, width, height);
    } else {
        return apply(im, width, height);
    }
}

void Resample::computeWeights(int oldSize, int newSize, vector<vector<pair<int, float> > > &matrix) {
    assert(newSize > 0, "Can only resample to positive sizes");

    float filterWidth = max(1.0f, (float)oldSize / newSize);

    matrix.resize(newSize);

    for (int x = 0; x < newSize; x++) {
        // This x in the output corresponds to which x in the input?
        float inX = (x + 0.5f) / newSize * oldSize - 0.5f;

        // Now compute a filter surrounding said x in the input
        int minX = ceilf(inX - filterWidth*3);
        int maxX = floorf(inX + filterWidth*3);
        minX = clamp(minX, 0, oldSize-1);
        maxX = clamp(maxX, 0, oldSize-1);

        assert(minX < maxX, "Wha?");

        // Compute a row of the sparse matrix
        matrix[x].resize(maxX - minX + 1);
        float totalWeight = 0;
        for (int i = minX; i <= maxX; i++) {
            float delta = i - inX;
            float w = lanczos_3(delta/filterWidth);
            matrix[x][i - minX] = make_pair(i, w);
            totalWeight += w;
        }
        for (int i = 0; i <= maxX - minX; i++) {
            matrix[x][i].second /= totalWeight;
        }
    }
}

Image Resample::resampleX(Image im, int width) {
    vector<vector<pair<int, float> > > matrix;
    computeWeights(im.width, width, matrix);

    Image out(width, im.height, im.frames, im.channels);

    for (int c = 0; c < out.channels; c++) {
        for (int t = 0; t < out.frames; t++) {
            for (int y = 0; y < out.height; y++) {
                for (int x = 0; x < out.width; x++) {
                    float val = 0;
                    for (size_t i = 0; i < matrix[x].size(); i++) {
                        val += matrix[x][i].second * im(matrix[x][i].first, y, t, c);
                    }
                    out(x, y, t, c) = val;
                }
            }
        }
    }

    return out;
}

Image Resample::resampleY(Image im, int height) {
    vector<vector<pair<int, float> > > matrix;
    computeWeights(im.height, height, matrix);

    Image out(im.width, height, im.frames, im.channels);

    for (int c = 0; c < out.channels; c++) {
        for (int t = 0; t < out.frames; t++) {
            for (int y = 0; y < out.height; y++) {
                for (int x = 0; x < out.width; x++) {
                    float val = 0;
                    for (size_t i = 0; i < matrix[y].size(); i++) {
                        val += matrix[y][i].second * im(x, matrix[y][i].first, t, c);
                    }
                    out(x, y, t, c) = val;
                }
            }
        }
    }

    return out;
}

Image Resample::resampleT(Image im, int frames) {
    vector<vector<pair<int, float> > > matrix;
    computeWeights(im.frames, frames, matrix);

    Image out(im.width, im.height, frames, im.channels);

    for (int c = 0; c < out.channels; c++) {
        for (int t = 0; t < out.frames; t++) {
            for (int y = 0; y < out.height; y++) {
                for (int x = 0; x < out.width; x++) {
                    float val = 0;
                    for (size_t i = 0; i < matrix[t].size(); i++) {
                        val += matrix[t][i].second * im(x, y, matrix[t][i].first, c);
                    }
                    out(x, y, t, c) = val;
                }
            }
        }
    }

    return out;
}



void Interleave::help() {
    pprintf("-interleave divides the image into n equally sized volumes and interleaves"
            " them. When given two arguments it operates on columns and rows. When"
            " given three arguments, it operates on columns, rows, and frames.\n\n"
            "Usage: ImageStack -load deck.jpg -interleave 2 2 -save shuffled.jpg\n\n");
}

bool Interleave::test() {
    Image a(120, 30, 20, 2);
    Noise::apply(a, -1, 1);
    Image b = a.copy();
    Interleave::apply(b, 2, 3, 4);
    for (int i = 0; i < 10; i++) {
        int tx = randomInt(0, 1);
        int ty = randomInt(0, 2);
        int tt = randomInt(0, 3);
        int x = randomInt(0, 59);
        int y = randomInt(0, 9);
        int t = randomInt(0, 4);
        for (int c = 0; c < a.channels; c++) {
            if (a(tx*60+x, ty*10+y, tt*5+t, c) !=
                b(x*2+tx, y*3+ty, t*4+tt, c)) {
                return false;
            }
        }
    }
    return true;
}

void Interleave::parse(vector<string> args) {
    if (args.size() == 2) {
        apply(stack(0), readInt(args[0]), readInt(args[1]));
    } else if (args.size() == 3) {
        apply(stack(0), readInt(args[0]), readInt(args[1]), readInt(args[2]));
    } else {
        panic("-interleave takes two or three arguments\n");
    }
}

void Interleave::apply(Image im, int rx, int ry, int rt) {
    assert(rt >= 1 && rx >= 1 && ry >= 1, "arguments to interleave must be strictly positive integers\n");

    // interleave in t
    if (rt != 1) {
        vector<float> tmp(im.frames);
        for (int c = 0; c < im.channels; c++) {
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width; x++) {
                    // copy out this chunk
                    for (int t = 0; t < im.frames; t++) {
                        tmp[t] = im(x, y, t, c);
                    }

                    // paste this chunk back in in bizarro order
                    int oldT = 0;
                    for (int t = 0; t < im.frames; t++) {
                        im(x, y, oldT, c) = tmp[t];
                        oldT += rt;
                        if (oldT >= im.frames) oldT = (oldT % rt) + 1;
                    }
                }
            }
        }
    }

    // interleave in x
    if (rx != 1) {
        vector<float> tmp(im.width);
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int y = 0; y < im.height; y++) {
                    for (int x = 0; x < im.width; x++) {
                        tmp[x] = im(x, y, t, c);
                    }

                    // paste this chunk back in in bizarro order
                    int oldX = 0;
                    for (int x = 0; x < im.width; x++) {
                        im(oldX, y, t, c) = tmp[x];
                        oldX += rx;
                        if (oldX >= im.width) oldX = (oldX % rx) + 1;
                    }
                }
            }
        }
    }

    // interleave in y
    if (ry != 1) {
        vector<float> tmp(im.height);
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int x = 0; x < im.width; x++) {
                    // copy out this chunk
                    for (int y = 0; y < im.height; y++) {
                        tmp[y] = im(x, y, t, c);
                    }

                    // paste this chunk back in in bizarro order
                    int oldY = 0;
                    for (int y = 0; y < im.height; y++) {
                        im(x, oldY, t, c) = tmp[y];
                        oldY += ry;
                        if (oldY >= im.height) oldY = (oldY % ry) + 1;
                    }
                }
            }
        }
    }
}

void Deinterleave::help() {
    pprintf("-deinterleave collects every nth frame, column, and/or row of the image"
            " and tiles the resulting collections. When given two arguments it"
            " operates on columns and rows. When given three arguments, it operates"
            " on all columns, rows, and frames.\n\n"
            "Usage: ImageStack -load lf.exr -deinterleave 16 16 -save lftranspose.exr\n\n");
}

bool Deinterleave::test() {
    Image a(120, 30, 20, 2);
    Noise::apply(a, -1, 1);
    Image b = a.copy();
    Deinterleave::apply(b, 2, 3, 4);
    for (int i = 0; i < 10; i++) {
        int tx = randomInt(0, 1);
        int ty = randomInt(0, 2);
        int tt = randomInt(0, 3);
        int x = randomInt(0, 59);
        int y = randomInt(0, 9);
        int t = randomInt(0, 4);
        for (int c = 0; c < a.channels; c++) {
            if (b(tx*60+x, ty*10+y, tt*5+t, c) !=
                a(x*2+tx, y*3+ty, t*4+tt, c)) {
                return false;
            }
        }
    }
    return true;
}

void Deinterleave::parse(vector<string> args) {
    if (args.size() == 2) {
        apply(stack(0), readInt(args[0]), readInt(args[1]));
    } else if (args.size() == 3) {
        apply(stack(0), readInt(args[0]), readInt(args[1]), readInt(args[2]));
    } else {
        panic("-deinterleave takes two or three arguments\n");
    }
}

void Deinterleave::apply(Image im, int rx, int ry, int rt) {
    assert(rt >= 1 && rx >= 1 && ry >= 1, "arguments to deinterleave must be strictly positive integers\n");

    // interleave in t
    if (rt != 1) {
        vector<float> tmp(im.frames);
        for (int c = 0; c < im.channels; c++) {
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width; x++) {
                    // copy out this chunk
                    for (int t = 0; t < im.frames; t++) {
                        tmp[t] = im(x, y, t, c);
                    }

                    // paste this chunk back in in bizarro order
                    int oldT = 0;
                    for (int t = 0; t < im.frames; t++) {
                        im(x, y, t, c) = tmp[oldT];
                        oldT += rt;
                        if (oldT >= im.frames) oldT = (oldT % rt) + 1;
                    }
                }
            }
        }
    }

    // interleave in x
    if (rx != 1) {
        vector<float> tmp(im.width);
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int y = 0; y < im.height; y++) {
                    for (int x = 0; x < im.width; x++) {
                        tmp[x] = im(x, y, t, c);
                    }

                    // paste this chunk back in in bizarro order
                    int oldX = 0;
                    for (int x = 0; x < im.width; x++) {
                        im(x, y, t, c) = tmp[oldX];
                        oldX += rx;
                        if (oldX >= im.width) oldX = (oldX % rx) + 1;
                    }
                }
            }
        }
    }

    // interleave in y
    if (ry != 1) {
        vector<float> tmp(im.height);
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int x = 0; x < im.width; x++) {
                    // copy out this chunk
                    for (int y = 0; y < im.height; y++) {
                        tmp[y] = im(x, y, t, c);
                    }

                    // paste this chunk back in in bizarro order
                    int oldY = 0;
                    for (int y = 0; y < im.height; y++) {
                        im(x, y, t, c) = tmp[oldY];
                        oldY += ry;
                        if (oldY >= im.height) oldY = (oldY % ry) + 1;
                    }
                }
            }
        }
    }

}


void Rotate::help() {
    printf("\n-rotate takes a number of degrees, and rotates every frame of the current image\n"
           "clockwise by that angle. The rotation preserves the image size, filling empty\n"
           " areas with zeros, and throwing away data which will not fit in the bounds.\n\n"
           "Usage: ImageStack -load a.tga -rotate 45 -save b.tga\n\n");
}

bool Rotate::test() {
    Image a(12, 13, 3, 2);
    Noise::apply(a, 4, 5);
    a = Resample::apply(a, 48, 52, 3);
    Image b = Rotate::apply(a, 70);
    b = Rotate::apply(b, 20);
    b = Rotate::apply(b, 80);
    b = Rotate::apply(b, 10);
    for (int i = 0; i < 100; i++) {
        int x = randomInt(10, a.width-11);
        int y = randomInt(10, a.height-11);
        int t = randomInt(0, a.frames-1);
        for (int c = 0; c < a.channels; c++) {
            if (!nearlyEqual(a(x, y, t, c),
                             b(a.width-1-x, a.height-1-y, t, c))) {
                return false;
            }
        }
    }
    return true;
}


void Rotate::parse(vector<string> args) {
    assert(args.size() == 1, "-rotate takes one argument\n");
    Image im = apply(stack(0), readFloat(args[0]));
    pop();
    push(im);
}


Image Rotate::apply(Image im, float degrees) {

    // figure out the rotation matrix
    float radians = degrees * M_PI / 180;
    float cosine = cos(radians);
    float sine = sin(radians);
    // locate the origin
    float xorigin = (im.width-1) * 0.5;
    float yorigin = (im.height-1) * 0.5;

    vector<float> matrix(6);
    matrix[0] = cosine; matrix[1] = sine; matrix[2] = xorigin - (cosine * xorigin + sine * yorigin);
    matrix[3] = -sine; matrix[4] = cosine; matrix[5] = yorigin - (-sine * xorigin + cosine * yorigin);
    return AffineWarp::apply(im, matrix);
}


void AffineWarp::help() {
    printf("\n-affinewarp takes a 2x3 matrix in row major order, and performs that affine warp\n"
           "on the image.\n\n"
           "Usage: ImageStack -load a.jpg -affinewarp 0.9 0.1 0 0.1 0.9 0 -save out.jpg\n\n");
}

bool AffineWarp::test() {
    // already tested by rotate
    return true;
}

void AffineWarp::parse(vector<string> args) {
    assert(args.size() == 6, "-affinewarp takes six arguments\n");
    vector<float> matrix(6);
    for (int i = 0; i < 6; i++) { matrix[i] = readFloat(args[i]); }
    Image im = apply(stack(0), matrix);
    pop();
    push(im);
}

Image AffineWarp::apply(Image im, vector<float> matrix) {

    assert(matrix.size() == 6, "An affine warp requires a vector with 6 entries\n");
    return apply(im, &matrix[0]);
}

Image AffineWarp::apply(Image im, float *matrix) {
    Image out(im.width, im.height, im.frames, im.channels);

    vector<float> sample(im.channels);
    for (int t = 0; t < im.frames; t++) {
        for (int y = 0; y < im.height; y++) {
            for (int x = 0; x < im.width; x++) {
                // figure out the sample location
                float fx = matrix[0] * x + matrix[1] * y + matrix[2];
                float fy = matrix[3] * x + matrix[4] * y + matrix[5];
                // don't sample outside the image
                if (fx < 0 || fx > im.width || fy < 0 || fy > im.height) {
                    for (int c = 0; c < im.channels; c++)
                        out(x, y, t, c) = 0;
                } else {
                    im.sample2D(fx, fy, t, sample);
                    for (int c = 0; c < im.channels; c++)
                        out(x, y, t, c) = sample[c];
                }
            }
        }
    }

    return out;
}

void Crop::help() {
    pprintf("-crop takes either zero, two, four, or six arguments. The first half"
            " of the arguments are either minimum t, minimum x and y, or all three"
            " in the order x, y, t. The second half of the arguments are"
            " correspondingly either number of frames, width and height, or all"
            " three in the order width, height, frames. You may crop outside the"
            " bounds of the original image. Values there are assumed to be black. If"
            " no argument are given, ImageStack guesses how much to crop by trimming"
            " rows and columns that are all the same color as the top left"
            " pixel.\n\n"
            "Usage: ImageStack -loadframes f*.tga -crop 10 1 -save frame10.tga\n"
            "       ImageStack -load a.tga -crop 100 100 200 200 -save cropped.tga\n"
            "       ImageStack -loadframes f*.tga -crop 100 100 10 200 200 1\n"
            "                  -save frame10cropped.tga\n\n");
}

bool Crop::test() {
    Image a(123, 234, 43, 2);
    Noise::apply(a, 0, 1);

    // within
    {
        Image b = Crop::apply(a, 2, 3, 4, 100, 200, 30);
        b -= a.region(2, 3, 4, 0, 100, 200, 30, 2);
        Stats sb(b);
        if (sb.mean() != 0 || sb.variance() != 0) return false;
    }

    // off the top left
    {
        Image b = Crop::apply(a, -5, -5, -5, 100, 200, 30);
        b.region(5, 5, 5, 0, 100-5, 200-5, 30-5, 2) -= a.region(0, 0, 0, 0, 100-5, 200-5, 30-5, 2);
        Stats sb(b);
        if (sb.mean() != 0 || sb.variance() != 0) return false;
    }

    // bottom right
    {
        Image b = Crop::apply(a, 50, 50, 10, 100, 200, 40);
        b.region(0, 0, 0, 0, 123-50, 234-50, 43-10, 2) -= a.region(50, 50, 10, 0, 123-50, 234-50, 43-10, 2);
        Stats sb(b);
        if (sb.mean() != 0 || sb.variance() != 0) return false;
    }

    return true;
}

void Crop::parse(vector<string> args) {

    Image im;

    if (args.size() == 0) {
        im = apply(stack(0));
    } else if (args.size() == 2) {
        im = apply(stack(0),
                   0, 0, readInt(args[0]),
                   stack(0).width, stack(0).height, readInt(args[1]));
    } else if (args.size() == 4) {
        im = apply(stack(0),
                   readInt(args[0]), readInt(args[1]),
                   readInt(args[2]), readInt(args[3]));
    } else if (args.size() == 6) {
        im = apply(stack(0),
                   readInt(args[0]), readInt(args[1]), readInt(args[2]),
                   readInt(args[3]), readInt(args[4]), readInt(args[5]));
    } else {
        panic("-crop takes two, four, or six arguments.\n");
    }

    pop();
    push(im);
}

Image Crop::apply(Image im) {
    int minX, maxX, minY, maxY, minT, maxT;

    // calculate minX
    for (minX = 0; minX < im.width; minX++) {
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int y = 0; y < im.height; y++) {
                    if (im(minX, y, t, c) != im(0, 0, 0, c)) { goto minXdone; }
                }
            }
        }
    }
minXdone:

    // calculate maxX
    for (maxX = im.width-1; maxX >= 0; maxX--) {
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int y = 0; y < im.height; y++) {
                    if (im(maxX, y, t, c) != im(0, 0, 0, c)) { goto maxXdone; }
                }
            }
        }
    }
maxXdone:

    // calculate minY
    for (minY = 0; minY < im.height; minY++) {
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int x = 0; x < im.width; x++) {
                    if (im(x, minY, t, c) != im(0, 0, 0, c)) { goto minYdone; }
                }
            }
        }
    }
minYdone:

    // calculate maxY
    for (maxY = im.height-1; maxY >= 0; maxY--) {
        for (int t = 0; t < im.frames; t++) {
            for (int x = 0; x < im.width; x++) {
                for (int c = 0; c < im.channels; c++) {
                    if (im(x, maxY, t, c) != im(0, 0, 0, c)) { goto maxYdone; }
                }
            }
        }
    }
maxYdone:

    // calculate minT
    for (minT = 0; minT < im.frames; minT++) {
        for (int y = 0; y < im.height; y++) {
            for (int x = 0; x < im.width; x++) {
                for (int c = 0; c < im.channels; c++) {
                    if (im(x, y, minT, c) != im(0, 0, 0, c)) { goto minTdone; }
                }
            }
        }
    }
minTdone:

    // calculate maxT
    for (maxT = im.frames-1; maxT >= 0; maxT--) {
        for (int y = 0; y < im.height; y++) {
            for (int x = 0; x < im.width; x++) {
                for (int c = 0; c < im.channels; c++) {
                    if (im(x, y, maxT, c) != im(0, 0, 0, c)) { goto maxTdone; }
                }
            }
        }
    }
maxTdone:

    int width = maxX - minX + 1;
    int height = maxY - minY + 1;
    int frames = maxT - minT + 1;

    assert(width >= 0 && height >= 0 && frames >= 0, "Can't auto crop a blank image\n");

    return apply(im, minX, minY, minT, width, height, frames);

}

Image Crop::apply(Image im, int minX, int minY, int width, int height) {
    return apply(im,
                 minX, minY, 0,
                 width, height, im.frames);
}


Image Crop::apply(Image im, int minX, int minY, int minT,
                  int width, int height, int frames) {
    Image out(width, height, frames, im.channels);

    for (int c = 0; c < im.channels; c++) {
        for (int t = max(0, -minT); t < min(frames, im.frames - minT); t++) {
            for (int y = max(0, -minY); y < min(height, im.height - minY); y++) {
                for (int x = max(0, -minX); x < min(width, im.width - minX); x++) {
                    out(x, y, t, c) = im(x + minX, y + minY, t + minT, c);
                }
            }
        }
    }

    return out;
}

void Flip::help() {
    printf("-flip takes 'x', 'y', or 't' as the argument and flips the current image along\n"
           "that dimension.\n\n"
           "Usage: ImageStack -load a.tga -flip x -save reversed.tga\n\n");
}

bool Flip::test() {
    Image a(123, 234, 43, 3);
    Noise::apply(a, -4, 3);
    Image fx = a.copy();
    Flip::apply(fx, 'x');
    Image fy = a.copy();
    Flip::apply(fy, 'y');
    Image ft = a.copy();
    Flip::apply(ft, 't');
    for (int i = 0; i < 10; i++) {
        int x = randomInt(0, a.width-1);
        int y = randomInt(0, a.height-1);
        int t = randomInt(0, a.frames-1);
        int c = randomInt(0, a.channels-1);
        if (a(x, y, t, c) != fx(a.width-1-x, y, t, c)) return false;
        if (a(x, y, t, c) != fy(x, a.height-1-y, t, c)) return false;
        if (a(x, y, t, c) != ft(x, y, a.frames-1-t, c)) return false;
    }
    return true;
}

void Flip::parse(vector<string> args) {
    assert(args.size() == 1, "-flip takes exactly one argument\n");
    char dimension = readChar(args[0]);
    apply(stack(0), dimension);
}

void Flip::apply(Image im, char dimension) {
    if (dimension == 't') {
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames/2; t++) {
                for (int y = 0; y < im.height; y++) {
                    for (int x = 0; x < im.width; x++) {
                        swap(im(x, y, t, c), im(x, y, im.frames-t-1, c));
                    }
                }
            }
        }
    } else if (dimension == 'y') {
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int y = 0; y < im.height/2; y++) {
                    for (int x = 0; x < im.width; x++) {
                        swap(im(x, y, t, c), im(x, im.height-1-y, t, c));
                    }
                }
            }
        }
    } else if (dimension == 'x') {
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int y = 0; y < im.height; y++) {
                    for (int x = 0; x < im.width/2; x++) {
                        swap(im(x, y, t, c), im(im.width-1-x, y, t, c));
                    }
                }
            }
        }
    } else {
        panic("-flip only understands dimensions 'x', 'y', and 't'\n");
    }
}


void Adjoin::help() {
    printf("\n-adjoin takes 'x', 'y', 't', or 'c' as the argument, and joins the top two\n"
           "images along that dimension. The images must match in the other dimensions.\n\n"
           "Usage: ImageStack -load a.tga -load b.tga -adjoin x -save ab.tga\n\n");
}

bool Adjoin::test() {
    Image a(32, 43, 23, 3);
    Noise::apply(a, 0, 1);

    Image fx(2, 43, 23, 3);
    Image fy(32, 2, 23, 3);
    Image ft(32, 43, 2, 3);
    Image fc(32, 43, 23, 2);
    Noise::apply(fx, 0, 1);
    Noise::apply(fy, 0, 1);
    Noise::apply(ft, 0, 1);
    Noise::apply(fc, 0, 1);
    Image afx = Adjoin::apply(fx, a, 'x');
    Image afy = Adjoin::apply(fy, a, 'y');
    Image aft = Adjoin::apply(ft, a, 't');
    Image afc = Adjoin::apply(fc, a, 'c');

    for (int i = 0; i < 10; i++) {
        int x = randomInt(0, a.width-1);
        int y = randomInt(0, a.height-1);
        int t = randomInt(0, a.frames-1);
        int c = randomInt(0, a.channels-1);
        int j = randomInt(0, 1);
        if (fx(j, y, t, c) != afx(j, y, t, c)) return false;
        if (fy(x, j, t, c) != afy(x, j, t, c)) return false;
        if (ft(x, y, j, c) != aft(x, y, j, c)) return false;
        if (fc(x, y, t, j) != afc(x, y, t, j)) return false;
        if (a(x, y, t, c) != afx(x+2, y, t, c)) return false;
        if (a(x, y, t, c) != afy(x, y+2, t, c)) return false;
        if (a(x, y, t, c) != aft(x, y, t+2, c)) return false;
        if (a(x, y, t, c) != afc(x, y, t, c+2)) return false;
    }
    return true;
}

void Adjoin::parse(vector<string> args) {
    assert(args.size() == 1, "-adjoin takes exactly one argument\n");
    char dimension = readChar(args[0]);
    Image im = apply(stack(1), stack(0), dimension);
    pop();
    pop();
    push(im);
}


Image Adjoin::apply(Image a, Image b, char dimension) {
    int newFrames = a.frames, newWidth = a.width, newHeight = a.height, newChannels = a.channels;
    int tOff = 0, xOff = 0, yOff = 0, cOff = 0;

    if (dimension == 't') {
        assert(a.width    == b.width &&
               a.height   == b.height &&
               a.channels == b.channels,
               "Cannot adjoin images that don't match in other dimensions\n");
        tOff = newFrames;
        newFrames += b.frames;
    } else if (dimension == 'y') {
        assert(a.width    == b.width &&
               a.frames   == b.frames &&
               a.channels == b.channels,
               "Cannot adjoin images that don't match in other dimensions\n");
        yOff = newHeight;
        newHeight += b.height;
    } else if (dimension == 'c') {
        assert(a.frames == b.frames &&
               a.height == b.height &&
               a.width  == b.width,
               "Cannot adjoin images that don't match in other dimensions\n");
        cOff = newChannels;
        newChannels += b.channels;
    } else if (dimension == 'x') {
        assert(a.frames == b.frames &&
               a.height == b.height &&
               a.channels  == b.channels,
               "Cannot adjoin images that don't match in other dimensions\n");
        xOff = newWidth;
        newWidth += b.width;
    } else {
        panic("-adjoin only understands dimensions 'x', 'y', and 't'\n");
    }

    Image out(newWidth, newHeight, newFrames, newChannels);
    // paste in the first image
    for (int c = 0; c < a.channels; c++) {
        for (int t = 0; t < a.frames; t++) {
            for (int y = 0; y < a.height; y++) {
                for (int x = 0; x < a.width; x++) {
                    out(x, y, t, c) = a(x, y, t, c);
                }
            }
        }
    }
    // paste in the second image
    for (int c = 0; c < b.channels; c++) {
        for (int t = 0; t < b.frames; t++) {
            for (int y = 0; y < b.height; y++) {
                for (int x = 0; x < b.width; x++) {
                    out(x + xOff, y + yOff, t + tOff, c + cOff) = b(x, y, t, c);
                }
            }
        }
    }

    return out;
}

void Transpose::help() {
    printf("-transpose takes two dimension of the form 'x', 'y', or 't' and transposes\n"
           "the current image over those dimensions. If given no arguments, it defaults\n"
           "to x and y.\n\n"
           "Usage: ImageStack -load a.tga -transpose x y -flip x -save rotated.tga\n\n");
}

void Transpose::parse(vector<string> args) {
    assert(args.size() == 0 || args.size() == 2, "-transpose takes either zero or two arguments\n");
    if (args.size() == 0) {
        Image im = apply(stack(0), 'x', 'y');
        pop();
        push(im);
    } else {
        char arg1 = readChar(args[0]);
        char arg2 = readChar(args[1]);
        Image im = apply(stack(0), arg1, arg2);
        pop();
        push(im);
    }
}

bool Transpose::test() {
    Image a(12, 23, 4, 2);
    Noise::apply(a, 12, 17);
    Image a1 = a;                        // 12 23 4 2
    a1 = Transpose::apply(a1, 'c', 'y'); // 12 2 4 23
    a1 = Transpose::apply(a1, 'x', 't'); // 4 2 12 23
    a1 = Transpose::apply(a1, 'c', 't'); // 4 2 23 12
    a1 = Transpose::apply(a1, 'y', 'x'); // 2 4 23 12

    Image a2 = a;                        // 12 23 4 2
    a2 = Transpose::apply(a2, 'x', 'c'); // 2 23 4 12
    a2 = Transpose::apply(a2, 'y', 't'); // 2 4 23 12

    a2 -= a1;
    Stats s(a2);
    if (s.mean() != 0 || s.variance() != 0) return false;

    return true;
}

Image Transpose::apply(Image im, char arg1, char arg2) {

    char dim1 = min(arg1, arg2);
    char dim2 = max(arg1, arg2);

    if (dim1 == 'c' && dim2 == 'y') {
        Image out(im.width, im.channels, im.frames, im.height);
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int y = 0; y < im.height; y++) {
                    for (int x = 0; x < im.width; x++) {
                        out(x, c, t, y) = im(x, y, t, c);
                    }
                }
            }
        }
        return out;
    } else if (dim1 == 'c' && dim2 == 't') {
        Image out(im.width, im.height, im.channels, im.frames);
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int y = 0; y < im.height; y++) {
                    for (int x = 0; x < im.width; x++) {
                        out(x, y, c, t) = im(x, y, t, c);
                    }
                }
            }
        }
        return out;
    } else if (dim1 == 'c' && dim2 == 'x') {
        Image out(im.channels, im.height, im.frames, im.width);
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int y = 0; y < im.height; y++) {
                    for (int x = 0; x < im.width; x++) {
                        out(c, y, t, x) = im(x, y, t, c);
                    }
                }
            }
        }
        return out;
    } else if (dim1 == 'x' && dim2 == 'y') {

        Image out(im.height, im.width, im.frames, im.channels);
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int y = 0; y < im.height; y++) {
                    for (int x = 0; x < im.width; x++) {
                        out(y, x, t, c) = im(x, y, t, c);
                    }
                }
            }
        }
        return out;

    } else if (dim1 == 't' && dim2 == 'x') {
        Image out(im.frames, im.height, im.width, im.channels);
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int y = 0; y < im.height; y++) {
                    for (int x = 0; x < im.width; x++) {
                        out(t, y, x, c) = im(x, y, t, c);
                    }
                }
            }
        }
        return out;
    } else if (dim1 == 't' && dim2 == 'y') {
        Image out(im.width, im.frames, im.height, im.channels);
        for (int c = 0; c < im.channels; c++) {
            for (int t = 0; t < im.frames; t++) {
                for (int y = 0; y < im.height; y++) {
                    for (int x = 0; x < im.width; x++) {
                        out(x, t, y, c) = im(x, y, t, c);
                    }
                }
            }
        }
        return out;
    } else {
        panic("-transpose only understands dimensions 'c', 'x', 'y', and 't'\n");
    }

    // keep compiler happy
    return Image();

}

void Translate::help() {
    pprintf("-translate moves the image data, leaving black borders. It takes two"
            " or three arguments. Two arguments are interpreted as a shift in x and"
            " a shift in y. Three arguments indicates a shift in x, y, and"
            " t. Negative values shift to the top left and positive ones to the"
            " bottom right. Fractional shifts are permitted; Lanczos sampling is"
            " used in this case.\n\n"
            "Usage: ImageStack -load in.jpg -translate -10 -10 -translate 20 20\n"
            "                  -translate -10 -10 -save in_border.jpg\n\n");
}

bool Translate::test() {
    Image a(2, 2, 2, 4);
    Noise::apply(a, 0, 1);
    a = Resample::apply(a, 120, 100, 20);
    Image b = a;
    b = Translate::apply(b, -17.5, 2.5, 0.5);
    b = Translate::apply(b, 10.0, -1.25, -0.25);
    b = Translate::apply(b, 7.5, -1.25, -0.25);

    b -= a;
    Stats s(b.region(40, 10, 10, 0,
                     40, 80, 1, 4));
    return (nearlyEqual(s.mean(), 0) && nearlyEqual(s.variance(), 0));
}

void Translate::parse(vector<string> args) {
    if (args.size() == 2) {
        Image im = apply(stack(0), readFloat(args[0]), readFloat(args[1]), 0);
        pop();
        push(im);
    } else if (args.size() == 3) {
        Image im = apply(stack(0), readFloat(args[0]), readFloat(args[1]), readFloat(args[2]));
        pop();
        push(im);
    } else {
        panic("-translate requires two or three arguments\n");
    }
}

Image Translate::apply(Image im, float xoff, float yoff, float toff) {
    Image current = im;
    Image out;

    // First do any non-integer translations
    if (xoff != floorf(xoff)) {
        out = applyX(im, xoff);
        current = out;
        xoff = 0;
    }

    if (yoff != floorf(yoff)) {
        out = applyY(current, yoff);
        current = out;
        yoff = 0;
    }

    if (toff != floorf(toff)) {
        out = applyT(current, toff);
        current = out;
        toff = 0;
    }

    // Now take care of the integer ones with a crop
    return Crop::apply(current, -xoff, -yoff, -toff, im.width, im.height, im.frames);
}

Image Translate::applyX(Image im, float xoff) {
    int ix = floorf(xoff);
    float fx = xoff - ix;
    // compute a 6-tap lanczos kernel
    float kernel[6];
    kernel[0] = lanczos_3(-3 + fx);
    kernel[1] = lanczos_3(-2 + fx);
    kernel[2] = lanczos_3(-1 + fx);
    kernel[3] = lanczos_3(0 + fx);
    kernel[4] = lanczos_3(1 + fx);
    kernel[5] = lanczos_3(2 + fx);
    double sum = 0;
    for (int i = 0; i < 6; i++) sum += kernel[i];
    for (int i = 0; i < 6; i++) kernel[i] /= sum;

    Image out(im.width, im.height, im.frames, im.channels);
    for (int t = 0; t < im.frames; t++) {
        for (int y = 0; y < im.height; y++) {
            for (int x = 0; x < im.width; x++) {
                for (int kx = -3; kx < 3; kx++) {
                    int imx = x - ix + kx;
                    if (imx < 0) { continue; }
                    if (imx >= im.width) { continue; }
                    float w = kernel[kx+3];
                    for (int c = 0; c < im.channels; c++) {
                        out(x, y, t, c) += w*im(imx, y, t, c);
                    }
                }
            }
        }
    }
    return out;
}

Image Translate::applyY(Image im, float yoff) {
    int iy = floorf(yoff);
    float fy = yoff - iy;
    // compute a 6-tap lanczos kernel
    float kernel[6];
    kernel[0] = lanczos_3(-3 + fy);
    kernel[1] = lanczos_3(-2 + fy);
    kernel[2] = lanczos_3(-1 + fy);
    kernel[3] = lanczos_3(0 + fy);
    kernel[4] = lanczos_3(1 + fy);
    kernel[5] = lanczos_3(2 + fy);
    double sum = 0;
    for (int i = 0; i < 6; i++) sum += kernel[i];
    for (int i = 0; i < 6; i++) kernel[i] /= sum;

    Image out(im.width, im.height, im.frames, im.channels);
    for (int t = 0; t < im.frames; t++) {
        for (int y = 0; y < im.height; y++) {
            for (int ky = -3; ky < 3; ky++) {
                int imy = y - iy + ky;
                if (imy < 0) { continue; }
                if (imy >= im.height) { continue; }
                float w = kernel[ky+3];
                for (int x = 0; x < im.width; x++) {
                    for (int c = 0; c < im.channels; c++) {
                        out(x, y, t, c) += w*im(x, imy, t, c);
                    }
                }
            }
        }
    }
    return out;
}

Image Translate::applyT(Image im, float toff) {
    int it = floorf(toff);
    float ft = toff - it;
    // compute a 6-tap lanczos kernel
    float kernel[6];
    kernel[0] = lanczos_3(-3 + ft);
    kernel[1] = lanczos_3(-2 + ft);
    kernel[2] = lanczos_3(-1 + ft);
    kernel[3] = lanczos_3(0 + ft);
    kernel[4] = lanczos_3(1 + ft);
    kernel[5] = lanczos_3(2 + ft);
    double sum = 0;
    for (int i = 0; i < 6; i++) sum += kernel[i];
    for (int i = 0; i < 6; i++) kernel[i] /= sum;

    Image out(im.width, im.height, im.frames, im.channels);
    for (int t = 0; t < im.frames; t++) {
        for (int kt = -3; kt < 3; kt++) {
            int imt = t - it + kt;
            if (imt < 0) { continue; }
            if (imt >= im.frames) { continue; }
            float w = kernel[kt+3];
            for (int y = 0; y < im.height; y++) {
                for (int x = 0; x < im.width; x++) {
                    for (int c = 0; c < im.channels; c++) {
                        out(x, y, t, c) += w*im(x, y, imt, c);
                    }
                }
            }
        }
    }
    return out;
}

void Paste::help() {
    printf("-paste places some of the second image in the stack inside the top image, at\n"
           "the specified location. -paste accepts two or three, six, or nine arguments.\n"
           "When given two or three arguments, it interprets these as x and y, or x, y,\n"
           "and t, and pastes the whole of the second image onto that location in the first\n"
           "image. If six or nine arguments are given, the latter four or six arguments\n"
           "specify what portion of the second image is copied. The middle two or three\n"
           "arguments specify the top left, and the last two or three arguments specify\n"
           "the size of the region to paste.\n\n"
           "The format is thus: -paste [desination origin] [source origin] [size]\n\n"
           "Usage: ImageStack -load a.jpg -push 820 820 1 3 -paste 10 10 -save border.jpg\n\n");
}

bool Paste::test() {
    Image a(100, 100, 20, 3);
    Noise::apply(a, 0, 12);
    Image orig = a.copy();
    Image b(20, 20, 10, 3);
    Noise::apply(b, 0, 12);
    Paste::apply(a, b, 10, 10, 5);
    b -= a.region(10, 10, 5, 0, 20, 20, 10, 3);
    Stats s(b);
    return (s.variance() == 0 && s.mean() == 0);
}

void Paste::parse(vector<string> args) {
    int xdst = 0, ydst = 0, tdst = 0;
    int xsrc = 0, ysrc = 0, tsrc = 0;
    int width = stack(1).width;
    int height = stack(1).height;
    int frames = stack(1).frames;

    if (args.size() == 2) {
        xdst = readInt(args[0]);
        ydst = readInt(args[1]);
    } else if (args.size() == 3) {
        xdst = readInt(args[0]);
        ydst = readInt(args[1]);
        tdst = readInt(args[2]);
    } else if (args.size() == 6) {
        xdst = readInt(args[0]);
        ydst = readInt(args[1]);
        xsrc = readInt(args[2]);
        ysrc = readInt(args[3]);
        width = readInt(args[4]);
        height = readInt(args[5]);
    } else if (args.size() == 9) {
        xdst = readInt(args[0]);
        ydst = readInt(args[1]);
        tdst = readInt(args[2]);
        xsrc = readInt(args[3]);
        ysrc = readInt(args[4]);
        tsrc = readInt(args[5]);
        width  = readInt(args[6]);
        height = readInt(args[7]);
        frames = readInt(args[8]);
    } else {
        panic("-paste requires two, three, six, or nine arguments\n");
    }

    apply(stack(0), stack(1),
          xdst, ydst, tdst,
          xsrc, ysrc, tsrc,
          width, height, frames);
    pull(1);
    pop();

}


void Paste::apply(Image into, Image from,
                  int xdst, int ydst,
                  int xsrc, int ysrc,
                  int width, int height) {
    apply(into, from,
          xdst, ydst, 0,
          xsrc, ysrc, 0,
          width, height, from.frames);
}

void Paste::apply(Image into, Image from,
                  int xdst, int ydst, int tdst) {
    apply(into, from,
          xdst, ydst, tdst,
          0, 0, 0,
          from.width, from.height, from.frames);
}

void Paste::apply(Image into, Image from,
                  int xdst, int ydst, int tdst,
                  int xsrc, int ysrc, int tsrc,
                  int width, int height, int frames) {
    assert(into.channels == from.channels,
           "Images must have the same number of channels\n");
    assert(tdst >= 0 &&
           ydst >= 0 &&
           xdst >= 0 &&
           tdst + frames <= into.frames &&
           ydst + height <= into.height &&
           xdst + width  <= into.width,
           "Cannot paste outside the target image\n");
    assert(tsrc >= 0 &&
           ysrc >= 0 &&
           xsrc >= 0 &&
           tsrc + frames <= from.frames &&
           ysrc + height <= from.height &&
           xsrc + width  <= from.width,
           "Cannot paste from outside the source image\n");
    for (int c = 0; c < into.channels; c++) {
        for (int t = 0; t < frames; t++) {
            for (int y = 0; y < height; y++) {
                for (int x = 0; x < width; x++) {
                    into(x + xdst, y + ydst, t + tdst, c) =
                        from(x + xsrc, y + ysrc, t + tsrc, c);
                }
            }
        }
    }
}

void Tile::help() {
    printf("\n-tile repeats the image along each dimension. It interprets two arguments as\n"
           "repetitions in x and y. Three arguments are interpreted as repetitions in x,\n"
           "y, and t.\n\n"
           "Usage: ImageStack -load a.tga -tile 2 2 -save b.tga\n\n");
}

bool Tile::test() {
    Image a(20, 20, 3, 2);
    Noise::apply(a, 0, 1);
    Image b = Tile::apply(a, 5, 4, 3);
    Stats sa(a), sb(b);
    return (nearlyEqual(sa.mean(), sb.mean()) &&
            nearlyEqual(sa.variance(), sb.variance()));
}

void Tile::parse(vector<string> args) {
    int tRepeat = 1, xRepeat = 1, yRepeat = 1;
    if (args.size() == 2) {
        xRepeat = readInt(args[0]);
        yRepeat = readInt(args[1]);
    } else if (args.size() == 3) {
        xRepeat = readInt(args[0]);
        yRepeat = readInt(args[1]);
        tRepeat = readInt(args[2]);
    } else {
        panic("-tile takes two or three arguments\n");
    }
    Image im = apply(stack(0), xRepeat, yRepeat, tRepeat);
    pop();
    push(im);
}

Image Tile::apply(Image im, int xRepeat, int yRepeat, int tRepeat) {

    Image out(im.width * xRepeat, im.height * yRepeat, im.frames * tRepeat, im.channels);

    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames * tRepeat; t++) {
            int imT = t % im.frames;
            for (int y = 0; y < im.height * yRepeat; y++) {
                int imY = y % im.height;
                for (int x = 0; x < im.width * xRepeat; x++) {
                    int imX = x % im.width;
                    out(x, y, t, c) = im(imX, imY, imT, c);
                }
            }
        }
    }

    return out;
}


void Subsample::help() {
    printf("\n-subsample subsamples the current image. Given two integer arguments, a and b,\n"
           "it selects one out of every a frames starting from frame b. Given four arguments,\n"
           "a, b, c, d, it selects one pixel out of every axb sized box, starting from pixel\n"
           "(c, d). Given six arguments, a, b, c, d, e, f, it selects one pixel from every\n"
           "axbxc volume, in the order width, height, frames starting at pixel (d, e, f).\n\n"
           "Usage: ImageStack -load in.jpg -subsample 2 2 0 0 -save smaller.jpg\n\n");
}

bool Subsample::test() {
    Image a(200, 200, 33, 2);
    Noise::apply(a, 0, 1);
    Image b = Subsample::apply(a, 2, 3, 1, 0, 1, 0);
    Stats sa(a), sb(b);
    return (nearlyEqual(sa.mean(), sb.mean()) &&
            nearlyEqual(sa.variance(), sb.variance()));
}

void Subsample::parse(vector<string> args) {
    if (args.size() == 2) {
        Image im = apply(stack(0), readInt(args[0]), readInt(args[1]));
        pop(); push(im);
    } else if (args.size() == 4) {
        Image im = apply(stack(0), readInt(args[0]), readInt(args[1]),
                         readInt(args[2]), readInt(args[3]));
        pop(); push(im);
    } else if (args.size() == 6) {
        Image im = apply(stack(0), readInt(args[0]), readInt(args[1]), readInt(args[2]),
                         readInt(args[3]), readInt(args[4]), readInt(args[5]));
        pop(); push(im);
    } else {
        panic("-subsample needs two, four, or six arguments\n");
    }
}

Image Subsample::apply(Image im, int boxFrames, int offsetT) {
    return apply(im, 1, 1, boxFrames, 0, 0, offsetT);
}

Image Subsample::apply(Image im, int boxWidth, int boxHeight,
                       int offsetX, int offsetY) {
    return apply(im, boxWidth, boxHeight, 1, offsetX, offsetY, 0);
}

Image Subsample::apply(Image im, int boxWidth, int boxHeight, int boxFrames,
                       int offsetX, int offsetY, int offsetT) {

    int newFrames = 0, newWidth = 0, newHeight = 0;
    for (int t = offsetT; t < im.frames; t += boxFrames) { newFrames++; }
    for (int x = offsetX; x < im.width;  x += boxWidth) { newWidth++; }
    for (int y = offsetY; y < im.height; y += boxHeight) { newHeight++; }

    Image out(newWidth, newHeight, newFrames, im.channels);

    for (int c = 0; c < im.channels; c++) {
        int outT = 0;
        for (int t = offsetT; t < im.frames; t += boxFrames) {
            int outY = 0;
            for (int y = offsetY; y < im.height; y += boxHeight) {
                int outX = 0;
                for (int x = offsetX; x < im.width; x += boxWidth) {
                    out(outX, outY, outT, c) = im(x, y, t, c);
                    outX++;
                }
                outY++;
            }
            outT++;
        }
    }

    return out;
}

void TileFrames::help() {
    printf("\n-tileframes takes a volume and lays down groups of frames in a grid, dividing\n"
           "the number of frames by the product of the arguments. It takes two arguments,\n"
           "the number of old frames across each new frame, and the number of frames down.\n"
           "each new frame. The first batch of frames will appear as the first row of the.\n"
           "first frame of the new volume.\n\n"
           "Usage: ImageStack -loadframes frame*.tif -tileframes 5 5 -saveframes sheet%%d.tif\n\n");
}

bool TileFrames::test() {
    Image a(20, 20, 20, 2);
    Noise::apply(a, 0, 1);
    Image b = TileFrames::apply(a, 5, 4);
    Stats sa(a), sb(b);
    return (nearlyEqual(sa.mean(), sb.mean()) &&
            nearlyEqual(sa.variance(), sb.variance()));
}

void TileFrames::parse(vector<string> args) {
    assert(args.size() == 2, "-tileframes takes two arguments\n");

    Image im = apply(stack(0), readInt(args[0]), readInt(args[1]));
    pop();
    push(im);
}

Image TileFrames::apply(Image im, int xTiles, int yTiles) {

    int newWidth = im.width * xTiles;
    int newHeight = im.height * yTiles;
    int newFrames = (int)(ceil((float)im.frames / (xTiles * yTiles)));

    Image out(newWidth, newHeight, newFrames, im.channels);

    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < newFrames; t++) {
            int outY = 0;
            for (int yt = 0; yt < yTiles; yt++) {
                for (int y = 0; y < im.height; y++) {
                    int outX = 0;
                    for (int xt = 0; xt < xTiles; xt++) {
                        int imT = (t * yTiles + yt) * xTiles + xt;
                        if (imT >= im.frames) { break; }
                        for (int x = 0; x < im.width; x++) {
                            out(outX, outY, t, c) = im(x, y, imT, c);
                            outX++;
                        }
                    }
                    outY++;
                }
            }
        }
    }

    return out;
}

void FrameTiles::help() {
    printf("\n-frametiles takes a volume where each frame is a grid and breaks each frame\n"
           "into multiple frames, one per grid element. The two arguments specify the grid\n"
           "size. This operation is the inverse of tileframes.\n\n"
           "Usage: ImageStack -loadframes sheet*.tif -frametiles 5 5 -saveframes frame%%d.tif\n\n");
}

bool FrameTiles::test() {
    Image a(20, 20, 20, 2);
    Noise::apply(a, -3, 3);
    Image b = TileFrames::apply(a, 5, 4);
    Image c = FrameTiles::apply(b, 5, 4);
    a -= c;
    Stats s(a);
    return (s.mean() == 0  && s.variance() == 0);
}


void FrameTiles::parse(vector<string> args) {
    assert(args.size() == 2, "-frametiles takes two arguments\n");

    Image im = apply(stack(0), readInt(args[0]), readInt(args[1]));
    pop();
    push(im);
}

Image FrameTiles::apply(Image im, int xTiles, int yTiles) {

    assert(im.width % xTiles == 0 &&
           im.height % yTiles == 0,
           "The image is not divisible by the given number of tiles\n");

    int newWidth = im.width / xTiles;
    int newHeight = im.height / yTiles;
    int newFrames = im.frames * xTiles * yTiles;

    Image out(newWidth, newHeight, newFrames, im.channels);

    for (int c = 0; c < im.channels; c++) {
        for (int t = 0; t < im.frames; t++) {
            int imY = 0;
            for (int yt = 0; yt < yTiles; yt++) {
                for (int y = 0; y < newHeight; y++) {
                    int imX = 0;
                    for (int xt = 0; xt < xTiles; xt++) {
                        int outT = (t * yTiles + yt) * xTiles + xt;
                        for (int x = 0; x < newWidth; x++) {
                            out(x, y, outT, c) = im(imX, imY, t, c);
                            imX++;
                        }
                    }
                    imY++;
                }
            }
        }
    }

    return out;
}


void Warp::help() {
    pprintf("-warp treats the top image of the stack as coordinates in the second"
            " image, and samples the second image accordingly. It takes no"
            " arguments. The number of channels in the top image is the"
            " dimensionality of the warp, and should be three or less.\n"
            "\n"
            "Usage: ImageStack -load in.jpg -push -evalchannels \"x+y\" \"y\" -warp -save out.jpg\n\n");
}

bool Warp::test() {
    Image a(100, 100, 1, 3);
    Noise::apply(a, 0, 1);
    // Do an affine warp in two ways
    vector<float> matrix(6);
    matrix[0] = 0.9; matrix[1] = 0.1; matrix[2] = 3;
    matrix[3] = -0.2; matrix[4] = 0.8; matrix[5] = 3;
    Image warped = AffineWarp::apply(a, matrix);
    Image warpField(100, 100, 1, 2);
    for (int y = 0; y < 100; y++) {
        for (int x = 0; x < 100; x++) {
            warpField(x, y, 0, 0) = matrix[0] * x + matrix[1] * y + matrix[2];
            warpField(x, y, 0, 1) = matrix[3] * x + matrix[4] * y + matrix[5];
        }
    }
    Image warped2 = Warp::apply(warpField, a);
    return nearlyEqual(warped, warped2);
}

void Warp::parse(vector<string> args) {
    assert(args.size() == 0, "warp takes no arguments\n");
    Image im = apply(stack(0), stack(1));
    pop();
    pop();
    push(im);
}

Image Warp::apply(Image coords, Image source) {

    Image out(coords.width, coords.height, coords.frames, source.channels);

    vector<float> sample(out.channels);
    if (coords.channels == 3) {
        for (int t = 0; t < coords.frames; t++) {
            for (int y = 0; y < coords.height; y++) {
                for (int x = 0; x < coords.width; x++) {
                    source.sample3D(coords(x, y, t, 0),
                                    coords(x, y, t, 1),
                                    coords(x, y, t, 2),
                                    sample);
                    for (int c = 0; c < out.channels; c++)
                        out(x, y, t, c) = sample[c];
                }
            }
        }
    } else if (coords.channels == 2) {
        for (int t = 0; t < coords.frames; t++) {
            for (int y = 0; y < coords.height; y++) {
                for (int x = 0; x < coords.width; x++) {
                    source.sample2D(coords(x, y, t, 0),
                                    coords(x, y, t, 1),
                                    t, sample);
                    for (int c = 0; c < out.channels; c++)
                        out(x, y, t, c) = sample[c];
                }
            }
        }
    } else {
        panic("index image must have two or three channels\n");
    }
    return out;
}



void Reshape::help() {
    printf("\n-reshape changes the way the memory of the current image is indexed. The four\n"
           "integer arguments specify a new width, height, frames, and channels.\n\n"
           "Usage: ImageStack -load movie.tmp -reshape width height*frames 1 channels\n"
           "                  -save filmstrip.tmp\n");
}

bool Reshape::test() {
    Image a(123, 23, 23, 2);
    Noise::apply(a, -4, 2);
    Image b = Reshape::apply(a, 2, 23, 123, 23);
    Stats sa(a), sb(b);
    return (nearlyEqual(sa.mean(), sb.mean()) &&
            nearlyEqual(sa.variance(), sb.variance()));
}

void Reshape::parse(vector<string> args) {
    assert(args.size() == 4, "-reshape takes four arguments\n");
    Image im = apply(stack(0),
                     readInt(args[0]), readInt(args[1]),
                     readInt(args[2]), readInt(args[3]));
    pop();
    push(im);

}

Image Reshape::apply(Image im, int x, int y, int t, int c) {
    assert(t *x *y *c == im.frames * im.width * im.height * im.channels,
           "New shape uses a different amount of memory that the old shape.\n");
    assert(im.dense(), "Input image is not densely packed in memory");
    Image out(x, y, t, c);
    memcpy(&out(0, 0, 0, 0), &im(0, 0, 0, 0), x*y*t*c*sizeof(float));
    return out;
}



#include "footer.h"
