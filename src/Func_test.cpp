#include "ImageStack.h"
#include "Func.h"

using namespace ImageStack;
using namespace ImageStack::Expr;

#define work(X) ((X+X*X*X)/(sqrt(X)+X*X))
//#define work(X) (X)

#define X_TILE_SIZE 256
#define Y_TILE_SIZE 32


void blur_fast(Image in, Image out) {
    __m256 one_third = _mm256_set1_ps(1.0f/3);

    for (int c = 0; c < in.channels; c++) {
        for (int t = 0; t < in.frames; t++) {


#pragma omp parallel for            
            for (int yTile = 0; yTile < in.height; yTile += Y_TILE_SIZE) {
                __m256 v0, v1, v2, sum, avg;
                float tmp[(X_TILE_SIZE)*(Y_TILE_SIZE+2)];
                for (int xTile = 0; xTile < in.width; xTile += X_TILE_SIZE) {
                    float *tmpPtr = (float *)tmp;
                    for (int y = -1; y < Y_TILE_SIZE+1; y++) {
                        const float *inPtr = &(in(xTile, yTile+y, t, c));
                        for (int x = 0; x < X_TILE_SIZE; x += 8) {          
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
                    for (int y = 0; y < Y_TILE_SIZE; y++) {
                        float *outPtr = &(out(xTile, yTile+y, t, c));
                        for (int x = 0; x < X_TILE_SIZE; x += 8) {
                            v0 = _mm256_loadu_ps(tmpPtr+(2*X_TILE_SIZE));
                            v1 = _mm256_loadu_ps(tmpPtr+X_TILE_SIZE);
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


Func blur_halide(Func in) {
    X x; Y y; C c;
    Func blurx = (in(x-1, y, c) + in(x, y, c) + in(x+1, y, c))/3;
    Func blury = (blurx(x, y-1, c) + blurx(x, y, c) + blurx(x, y+1, c))/3;
    return blury;
}



Func blur_halide2(Func in) {
    Func blurx = (shiftX(in, -1) + in + shiftX(in, +1))/3;
    Func blury = (shiftY(blurx, -1) + blurx + shiftY(blurx, +1))/3;
    return blury;
}



int main(int argc, char **argv) {
    start();

    Image input = Load::apply(argv[1]);
    input = input.selectColumns(0, ((input.width-2)/X_TILE_SIZE)*X_TILE_SIZE+2);
    input = input.selectRows(0, ((input.height-2)/Y_TILE_SIZE)*Y_TILE_SIZE+2);

    printf("Using %d x %d of the input\n", input.width, input.height);

    const int iterations = 20;

    try {

        Image noise(128, 128, 128, 1);
        Noise::apply(noise, 0, 1);
        Image testY = interleaveY(noise, 0);
        Save::apply(testY, "interleaveY.tmp");

        Image output(input.width, input.height, input.frames, input.channels);
        double t;

        Func f = input+1;
        output = f;
        
        output.set(0);
        t = 1e10;
        for (int i = 0; i < iterations; i++) {
            double t1 = currentTime();
            output.set(blur_halide(zeroBoundary(input)));
            t = std::min(t, currentTime() - t1);
        }
        printf("%f\n", t);
        Save::apply(output, "output1.tmp");
        
        output.set(0);        
        t = 1e10;
        for (int i = 0; i < iterations; i++) {
            double t1 = currentTime();
            output.set(blur_halide2(zeroBoundary(input)));
            t = std::min(t, currentTime() - t1);
        }
        printf("%f\n", t);
        Save::apply(output, "output2.tmp");

        
        output.set(0);
        t = 1e10;
        for (int i = 0; i < iterations; i++) {
            double t1 = currentTime();
            auto zb = zeroBoundary(input);
            Func blurX = (shiftX(zb, -1) + zb + shiftX(zb, 1))/3;           
            output.set((shiftY(blurX, -1) + blurX + shiftY(blurX, 1))/3);
            t = std::min(t, currentTime() - t1);
        }
        printf("%f\n", t);
        Save::apply(output, "output3.tmp");

        output.set(0);
        t = 1e10;
        for (int i = 0; i < iterations; i++) {
            double t1 = currentTime();        
            blur_fast(input.region(1, 1, 0, 0, input.width-2, input.height-2, input.frames, input.channels),
                      output.region(1, 1, 0, 0, input.width-2, input.height-2, input.frames, input.channels));
            t = std::min(t, currentTime() - t1);
        }
        printf("%f\n", t);
        Save::apply(output, "output4.tmp");
        

    } catch (Exception &e) {
        printf("Failure: %s\n", e.message);
    }

    return 0;
}





