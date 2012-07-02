#ifndef IMAGESTACK_GKDTREE_H
#define IMAGESTACK_GKDTREE_H
#include "header.h"

#include <limits>
#include <math.h>
#include <stdlib.h>
#include <string.h>

float rand_float() {
    return rand()/(RAND_MAX+1.0f);
}

inline float gCDF(float x) {
    x *= 0.81649658092772592f;
    if (x < -2) {
        return 0;
    }
    if (x < -1) {
        x += 2;
        x *= x;
        x *= x;
        return x;
    }
    if (x < 0) {
        return 12 + x*(16 - x*x*(8 + 3*x));
    }
    if (x < 1) {
        return 12 + x*(16 - x*x*(8 - 3*x));
    }
    if (x < 2) {
        x = x-2;
        x *= x;
        x *= x;
        return -x + 24;
    }
    return 24;
}

class GKDTree {
public:
    /*
      static Image filter(Image im, Image ref) {
          Image out(im.width, im.height, im.frames, im.channels);

          vector<float *> points(ref.width*ref.height*ref.frames);
          int i = 0;
          for (int t = 0; t < ref.frames; t++) {
              for (int x = 0; x < ref.width; x++) {
                  for (int y = 0; y < ref.height; y++) {
                      points[i++] = ref(x, y, t);
                  }
              }
          }

          GKDTree tree(ref.channels, &points[0], points.size(), 2*0.707107);
          tree.finalize();

          vector<int> indices(64);
          vector<float> weights(64);

          Image leafValues(tree.getLeaves(), 1, 1, im.channels+1);

          float *imPtr = im(0, 0, 0);
          float *refPtr = ref(0, 0, 0);
          for (int t = 0; t < im.frames; t++) {
              for (int y = 0; y < im.height; y++) {
                  for (int x = 0; x < im.width; x++) {
                      int results = tree.gaussianLookup(refPtr,
                                                        &indices[0],
                                                        &weights[0],
                                                        4);
                      for (int i = 0; i < results; i++) {
                          float w = weights[i];
                          float *vPtr = leafValues(indices[i], 0);
                          for (int c = 0; c < im.channels; c++) {
                              vPtr[c] += imPtr[c]*w;
                          }
                          vPtr[im.channels] += w;
                      }
                      refPtr += ref.channels;
                      imPtr += im.channels;
                  }
              }
          }

          float *outPtr = out(0, 0, 0);
          float *slicePtr = ref(0, 0, 0);

          for (int t = 0; t < out.frames; t++) {
              for (int y = 0; y < out.height; y++) {
                  for (int x = 0; x < out.width; x++) {
                      int results = tree.gaussianLookup(&slicePtr[0],
                                                        &indices[0],
                                                        &weights[0],
                                                        64);
                      float outW = 0;

                      for (int i = 0; i < results; i++) {
                          float w = weights[i];
                          float *vPtr = leafValues(indices[i], 0);
                          for (int c = 0; c < out.channels; c++) {
                              outPtr[c] += vPtr[c]*w;
                          }
                          outW += w*vPtr[out.channels];
                      }

                      if (outW < 0.00000001) {
                          for (int c = 0; c < out.channels; c++) {
                              outPtr[c] = im(x, y, t)[c];
                          }
                      } else {
                          float invOutW = 1.0f/outW;
                          for (int c = 0; c < out.channels; c++) {
                              outPtr[c] *= invOutW;
                          }
                      }

                      slicePtr += ref.channels;
                      outPtr += out.channels;
                  }
              }
          }

          return out;
      }
    */

    // Build a gkdtree using the supplied array of points to control
    // the sampling.  sizeBound specifies the maximum allowable side
    // length of a kdtree leaf.  At least one point from data lies in
    // any given leaf.

    GKDTree(int dims, float **data, int nData, float sBound) :
        dimensions(dims), sizeBound(sBound), leaves(0) {

        root = build(data, nData);
    }

    ~GKDTree() {
        delete root;
    }

    void finalize() {
        float *kdtreeMins = new float[dimensions];
        float *kdtreeMaxs = new float[dimensions];

        for (int i = 0; i < dimensions; i++) {
            kdtreeMins[i] = -INF;
            kdtreeMaxs[i] = +INF;
        }

        root->computeBounds(kdtreeMins, kdtreeMaxs);

        delete[] kdtreeMins;
        delete[] kdtreeMaxs;
    }

    int getLeaves() {
        return leaves;
    }

    // Compute a gaussian spread of kdtree leaves around the given
    // point. This is the general case sampling strategy.
    int gaussianLookup(float *value, int *ids, float *weights, int nSamples) {
        return root->gaussianLookup(value, &ids, &weights, nSamples, 1);
    }

private:

    class Node {
    public:
        virtual ~Node() {}

        // Returns a list of samples from the kdtree distributed
        // around value with std-dev sigma in all dimensions. Some
        // samples may be repeated. Returns how many entries in the
        // ids and weights arrays were used.
        virtual int gaussianLookup(float *value, int **ids, float **weights, int nSamples, float p) = 0;

        // special case optimization of the above where nsamples = 1
        virtual int singleGaussianLookup(float *value, int **ids, float **weights, float p) = 0;

        virtual void computeBounds(float *mins, float *maxs) = 0;

    };

    class Split : public Node {
    public:
        virtual ~Split() {
            delete left;
            delete right;
        }


        // for a given gaussian and a given value, the probability of splitting left at this node
        inline float pLeft(float value) {
            // Coarsely approximate the cumulative normal distribution
            float val = gCDF(cut_val - value);
            float minBound = gCDF(min_val - value);
            float maxBound = gCDF(max_val - value);
            return (val - minBound) / (maxBound - minBound);
        }

        int gaussianLookup(float *value, int **ids, float **weights, int nSamples, float p) {
            // Calculate how much of a gaussian ball of radius sigma,
            // that has been trimmed by all the cuts so far, lies on
            // each side of the split

            // compute the probability of a sample splitting left
            float val = pLeft(value[cut_dim]);

            // Send some samples to the left of the split
            int leftSamples = (int)(val*nSamples);

            // Send some samples to the right of the split
            int rightSamples = (int)((1-val)*nSamples);

            // There's probably one sample left over by the rounding
            if (leftSamples + rightSamples != nSamples) {
                float fval = val*nSamples - leftSamples;
                // if val is high we send it left, if val is low we send it right
                if (rand_float() < fval) {
                    leftSamples++;
                } else {
                    rightSamples++;
                }
            }

            int samplesFound = 0;
            // Get the left samples
            if (leftSamples > 0) {
                if (leftSamples > 1) {
                    samplesFound += left->gaussianLookup(value, ids, weights, leftSamples, p*val);
                } else {
                    samplesFound += left->singleGaussianLookup(value, ids, weights, p*val);
                }
            }

            // Get the right samples
            if (rightSamples > 0) {
                if (rightSamples > 1) {
                    samplesFound += right->gaussianLookup(value, ids, weights, rightSamples, p*(1-val));
                } else {
                    samplesFound += right->singleGaussianLookup(value, ids, weights, p*(1-val));
                }
            }

            return samplesFound;
        }

        // a special case optimization of the above for when nSamples is 1
        int singleGaussianLookup(float *value, int **ids, float **weights, float p) {
            float val = pLeft(value[cut_dim]);
            if (rand_float() < val) {
                return left->singleGaussianLookup(value, ids, weights, p*val);
            } else {
                return right->singleGaussianLookup(value, ids, weights, p*(1-val));
            }
        }

        void computeBounds(float *mins, float *maxs) {
            min_val = mins[cut_dim];
            max_val = maxs[cut_dim];

            maxs[cut_dim] = cut_val;
            left->computeBounds(mins, maxs);
            maxs[cut_dim] = max_val;

            mins[cut_dim] = cut_val;
            right->computeBounds(mins, maxs);
            mins[cut_dim] = min_val;
        }

        int cut_dim;
        float cut_val, min_val, max_val;
        Node *left, *right;
    };

    class Leaf : public Node {

    public:
        Leaf(int id_, float **data, int nData, int dimensions_)
            : id(id_), dimensions(dimensions_) {
            position = new float[dimensions];
            for (int i = 0; i < dimensions; i++) {
                position[i] = 0;
                for (int j = 0; j < nData; j++) {
                    position[i] += data[j][i];
                }
                position[i] /= nData;
            }
        }

        ~Leaf() {
            delete[] position;
        }

        int gaussianLookup(float *query, int **ids, float **weights, int nSamples, float p) {
            // p is the probability with which one sample arrived here
            // calculate the correct probability, q

            float q = 0;
            for (int i = 0; i < dimensions; i++) {
                float diff = query[i] - position[i];
                diff *= diff;
                q += diff;
            }

            // Gaussian of variance 1/2
            q = expf(-q);

            *(*ids)++ = id;
            *(*weights)++ = nSamples * q / p;

            return 1;
        }

        int singleGaussianLookup(float *query, int **ids, float **weights, float p) {
            return gaussianLookup(query, ids, weights, 1, p);
        }

        void computeBounds(float *mins, float *maxs) {
        }

    private:
        int id, dimensions;
        float *position;
    };

    Node *root;
    int dimensions;
    float sizeBound;
    int leaves;

    Node *build(float **data, int nData) {

        if (nData == 1) {
            return new Leaf(leaves++, data, nData, dimensions);
        } else {

            vector<float> mins(dimensions), maxs(dimensions);

            // calculate the data bounds in every dimension
            for (int i = 0; i < dimensions; i++) {
                mins[i] = maxs[i] = data[0][i];
            }
            for (int j = 1; j < nData; j ++) {
                for (int i = 0; i < dimensions; i++) {
                    if (data[j][i] < mins[i]) mins[i] = data[j][i];
                    if (data[j][i] > maxs[i]) maxs[i] = data[j][i];
                }
            }

            // find the longest dimension
            int longest = 0;
            for (int i = 1; i < dimensions; i++) {
                float delta = maxs[i] - mins[i];
                if (delta > maxs[longest] - mins[longest])
                    longest = i;
            }

            // if it's large enough, cut in that dimension
            if (maxs[longest] - mins[longest] > sizeBound) {
                Split *n = new Split;
                n->cut_dim = longest;
                n->cut_val = (maxs[longest] + mins[longest])/2;

                // these get computed later
                n->min_val = -INF;
                n->max_val = INF;

                // resort the input over the split
                int pivot = 0;
                for (int i = 0; i < nData; i++) {
                    // The next value is larger than the pivot
                    if (data[i][longest] >= n->cut_val) continue;

                    // We haven't seen anything larger than the pivot yet
                    if (i == pivot) {
                        pivot++;
                        continue;
                    }

                    // The current value is smaller than the pivot
                    float *tmp = data[i];
                    data[i] = data[pivot];
                    data[pivot] = tmp;
                    pivot++;
                }

                // Build the two subtrees
                n->left = build(data, pivot);

                n->right = build(data+pivot, nData-pivot);

                return n;
            } else {
                return new Leaf(leaves++, data, nData, dimensions);
            }
        }
    };


};


#include "footer.h"

#endif
