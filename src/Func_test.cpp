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


int main(int argc, char **argv) {
    start();

    try {
        Func f;
        X x; Y y;
        f = x;
        Image foo(128, 128, 1, 1);        
        foo.set(f(toInt(sqrt(toFloat(x+y))), x-y, 0, 0));
        Save::apply(foo, "foo.tmp");
    } catch (Exception &e) {
        printf("Failure: %s\n", e.message);
    }

    return 0;
}





