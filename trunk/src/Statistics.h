#ifndef IMAGESTACK_STATISTICS_H
#define IMAGESTACK_STATISTICS_H
#include "header.h"

class Dimensions : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
};

class Stats {
public:
    Stats(Image im);

#define BASIC if (!basicStatsComputed) computeBasicStats();
#define MOMENT if (!momentsComputed) computeMoments();

    inline double sum(int c)     {BASIC; return sums[c];}
    inline double sum()          {BASIC; return sum_;}
    inline double mean(int c)    {BASIC; return means[c];}
    inline double mean()         {BASIC; return mean_;}
    inline double minimum(int c) {BASIC; return mins[c];}
    inline double minimum()      {BASIC; return min_;}
    inline double maximum(int c) {BASIC; return maxs[c];}
    inline double maximum()      {BASIC; return max_;}
    inline int nans() {BASIC; return nans_;}
    inline int posinfs() {BASIC; return posinfs_;}
    inline int neginfs() {BASIC; return neginfs_;}
    inline double covariance(int c1, int c2) {MOMENT; return covarianceMatrix[c1 * channels + c2];}
    inline double variance(int c) {MOMENT; return variances[c];}
    inline double variance()      {MOMENT; return variance_;}
    inline double skew(int c)     {MOMENT; return skews[c];}
    inline double skew()          {MOMENT; return skew_;}
    inline double kurtosis(int c) {MOMENT; return kurtoses[c];}
    inline double kurtosis()      {MOMENT; return kurtosis_;}
    inline double barycenterX(int c) { MOMENT; return barycenters[c*2]; }
    inline double barycenterY(int c) { MOMENT; return barycenters[c*2+1]; }
    inline double spatialVarianceX(int c) { MOMENT; return spatialVariances[c*2]; }
    inline double spatialVarianceY(int c) { MOMENT; return spatialVariances[c*2+1]; }

#undef BASIC
#undef MOMENT

private:
    void computeBasicStats();
    bool basicStatsComputed;
    void computeMoments();
    bool momentsComputed;
    Image im_;

    int channels;
    vector<double> sums, means, variances, kurtoses, skews, mins, maxs;
    vector<double> barycenters, spatialVariances;
    vector<double> covarianceMatrix;
    double sum_, mean_, variance_, min_, max_, kurtosis_, skew_;
    int nans_, neginfs_, posinfs_;
};

class Statistics : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im);
};

class Noise : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, float minVal, float maxVal);
};

class Histogram : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, int buckets = 256, float minVal = 0, float maxVal = 1);
};


class Equalize : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, float lower, float upper);
};


class HistogramMatch : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, Image model);
};


class Shuffle : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im);
};

class KMeans : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, int clusters);
};

class Sort : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, char dimension);
};

class DimensionReduction : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, int newChannels);
};

class LocalMaxima : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);

    struct Maximum {
        float x, y, t, value;

        Maximum() : x(0), y(0), t(0), value(0) {}
        Maximum(float x_, float y_, float t_, float value_) :
            x(x_), y(y_), t(t_), value(value_) {}

        bool operator<(const Maximum &other) const {
            return (value < other.value);
        }

        bool operator>(const Maximum &other) const {
            return (value > other.value);
        }

        bool operator<=(const Maximum &other) const {
            return (value <= other.value);
        }

        bool operator>=(const Maximum &other) const {
            return (value >= other.value);
        }

    };
    static vector<Maximum> apply(Image im, bool xCheck, bool yCheck, bool tCheck, float threshold, float minDistance);
};

class Printf : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, string fmt, vector<float> args);
};

class FPrintf : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static void apply(Image im, string filename, string fmt, vector<float> args);
};

class PCA : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, int newChannels);
};

class PatchPCA : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, float sigma, int newChannels);
};

class PatchPCA3D : public Operation {
public:
    void help();
    bool test();
    void parse(vector<string> args);
    static Image apply(Image im, float sigma, int newChannels);
};

#include "footer.h"
#endif
