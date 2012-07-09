#include "ImageStack.h"
#include "Func.h"

using namespace ImageStack;
using namespace ImageStack::Lazy;

//#define work(X) ((X+X*X*X)/(sqrt(X)+X*X))
#define work(X) (X)

void blur_fast(Image in, Image out) {
    __m256 one_third = _mm256_set1_ps(1.0f/3);

    for (int c = 0; c < in.channels; c++) {
        for (int t = 0; t < in.frames; t++) {
#pragma omp parallel for
            
            for (int yTile = 0; yTile < in.height; yTile += 64) {
                __m256 v0, v1, v2, sum, avg;
                float tmp[(64)*(64+2)];
                for (int xTile = 0; xTile < in.width; xTile += 64) {
                    float *tmpPtr = (float *)tmp;
                    for (int y = -1; y < 64+1; y++) {
                        const float *inPtr = &(in(xTile, yTile+y, t, c));
                        for (int x = 0; x < 64; x += 8) {          
                            v0 = _mm256_loadu_ps(inPtr-1);
                            v1 = _mm256_loadu_ps(inPtr+1);
                            v2 = _mm256_loadu_ps(inPtr);
                            sum = _mm256_add_ps(_mm256_add_ps(v0, v1), v2);
                            avg = _mm256_mul_ps(sum, one_third);
                            _mm256_storeu_ps(tmpPtr, avg);
                            tmpPtr += 8;
                            inPtr += 8;
                        }
                    }
                    tmpPtr = (float *)tmp;
                    for (int y = 0; y < 64; y++) {
                        float *outPtr = &(out(xTile, yTile+y, t, c));
                        for (int x = 0; x < 64; x += 8) {
                            v0 = _mm256_loadu_ps(tmpPtr+(2*64));
                            v1 = _mm256_loadu_ps(tmpPtr+64);
                            v2 = _mm256_loadu_ps(tmpPtr);
                            tmpPtr += 8;
                            sum = _mm256_add_ps(_mm256_add_ps(v0, v1), v2);
                            avg = _mm256_mul_ps(sum, one_third);
                            _mm256_storeu_ps(outPtr, avg);
                            outPtr += 8;
                        }
                    }
                } 
            }  
        }
    }
}


int main(int argc, char **argv) {
    start();

    try {
        Image in = Load::apply("in.tmp");       
        Image out(in.width, in.height, in.frames, in.channels);

        printf("Inline...\n");
        out.set(0);
        double t1s = currentTime();
        for (int i = 0; i < 10; i++) {
            // blur the square root
            // inline
            auto tmp = zeroBoundary(work(in));
            auto bx = shiftX(tmp, -1) + tmp + shiftX(tmp, 1);
            auto by = shiftY(bx, -1) + bx + shiftY(bx, 1);
            out.set(by*(1.0f/9));
        }
        double t1e = currentTime();
        Save::apply(out, "out_inline.tmp");


        printf("Root...\n");
        out.set(0);
        double t2s = currentTime();
        for (int i = 0; i < 10; i++) {
            // root
            Image tmp = work(in);
            Image bx = shiftX(zeroBoundary(tmp), -1) + tmp + shiftX(zeroBoundary(tmp), 1);
            Image by = (shiftY(zeroBoundary(bx), -1) + bx + shiftY(zeroBoundary(bx), 1))/9;
            out = by;
        }
        double t2e = currentTime();
        Save::apply(out, "out_root.tmp");

        printf("Chunk...\n");        
        out.set(0);
        double t3s = currentTime();
        for (int i = 0; i < 10; i++) {
            // chunk
            Func tmp = zeroBoundary(work(in));
            Func bx = shiftX(tmp, -1) + tmp + shiftX(tmp, 1);
            Func by = (shiftY(bx, -1) + bx + shiftY(bx, 1))/9;
            out.set(by);
        }
        double t3e = currentTime();
        Save::apply(out, "out_chunk.tmp");

        printf("Manual...\n");
        out.set(0);
        double t4s = currentTime();
        for (int i = 0; i < 10; i++) {
            blur_fast(in.region(1, 1, 0, 0, in.width-2, in.height-2, in.frames, in.channels), 
                      out.region(1, 1, 0, 0, in.width-2, in.height-2, in.frames, in.channels));
        }
        double t4e = currentTime();
        Save::apply(out, "out_manual.tmp");


        printf("%f %f %f %f\n", t1e-t1s, t2e-t2s, t3e-t3s, t4e-t4s);
        
    } catch (Exception &e) {
        printf("Failure: %s\n", e.message);
    }

    return 0;
}





