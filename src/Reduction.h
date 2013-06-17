#ifndef IMAGESTACK_REDUCTION_H
#define IMAGESTACK_REDUCTION_H

#include "Expr.h"
namespace ImageStack {

namespace Reduce {

class RowSum {
public:
    void accept(int, float val) {
        accum_scalar += val;
    }

    void accept(int, Vec::type val) {
        accum_vec += val;
    }    

    float toScalar() {
        float result = accum_scalar;
        float vec[Vec::width];
        Vec::store(accum_vec, vec);
        for (int i = 0; i < Vec::width; i++) {
            result += vec[i];
        }
        return result;
    }

    RowSum() : accum_vec(Vec::zero()), accum_scalar(0) {
    }


    Vec::type accum_vec;
    float accum_scalar;

};

template<typename T>
double sum(const T expr_, const FloatExprType(T) *ptr = NULL) {

    FloatExprType(T) expr(expr_);
        
    // Compute the domain over which we can vectorize
    bool boundedVX = expr.boundedVecX();
    int minVX = expr.minVecX();
    int maxVX = expr.maxVecX();
        
    int width = expr.getSize(0);
    int height = expr.getSize(1);
    int frames = expr.getSize(2);
    int channels = expr.getSize(3);

    Expr::Region r = {0, 0, 0, 0, width, height, frames, channels};

    // Figure out what regions of what functions are required here
    expr.prepare(r, 0);
    expr.prepare(r, 1);
    expr.prepare(r, 2);

    vector<float> rowSums(height);
    double total = 0.0;

    for (int c = 0; c < channels; c++) {
        for (int t = 0; t < frames; t++) {            

            #ifdef _OPENMP
            #pragma omp parallel for
            #endif            
            for (int y = 0; y < height; y++) {
                RowSum rowSum;
                FloatExprType(T)::Iter iter = expr.scanline(0, y, t, c, width);
                Expr::evaluateInto(iter, rowSum, 0, width, boundedVX, minVX, maxVX);                
                rowSums[y] = rowSum.toScalar();
            }

            for (int y = 0; y < height; y++) {
                total += rowSums[y];
            }
            
        }               
    }

    expr.prepare(r, 3);

    return total;
}

};

}

#endif
