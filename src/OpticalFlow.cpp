
// Implementation of an optical flow algorithm
// High accuracy optical flow estimation based on a theory for warping
// Algorithm from T. Brox et al. ECCV 2004.
// Implementation from C. Liu. ( http://people.csail.mit.edu/celiu/OpticalFlow/ )
// Ported to ImageStack by Sung Hee Park
// This code is for research purposes only. Contact authors for commercial use.

#include "main.h"
#include "File.h"
#include "Geometry.h"
#include "OpticalFlow.h"
#include "Arithmetic.h"
#include "Calculus.h"
#include "Statistics.h"
#include "Filter.h"
#include "Convolve.h"
#include "Color.h"
#include "Paint.h"
#include "header.h"

class GaussianPyramid {

private:
    Image *ImPyramid;
    int nLevels;

public:
    GaussianPyramid(void) {
        ImPyramid = NULL;
    }
    ~GaussianPyramid(void) {
        if (ImPyramid != NULL) {
            delete [] ImPyramid;
        }
    }

    void ConstructPyramid(Image image, double ratio = 0.8, int minWidth = 10) {

        // the ratio cannot be arbitrary numbers
        if (ratio > 0.98 || ratio < 0.4) {
            ratio = 0.75;
        }
        // first decide how many levels
        nLevels=log((double)minWidth/image.width) / log(ratio);
        if (ImPyramid != NULL) {
            delete [] ImPyramid;
        }
        ImPyramid = new Image[nLevels];
        ImPyramid[0] = image.copy();
        double baseSigma = (1.0/ratio - 1.0);
        for (int i=1; i < nLevels; i++) {
            Image foo = image.copy();
            double sigma = baseSigma * i;
            double scale = pow(ratio, i);
            FastBlur::apply(foo, sigma, sigma, 0);
            ImPyramid[i] = Resample::apply(foo, foo.width*scale, foo.height*scale);
        }
    }

    inline int nlevels() const {return nLevels;};
    inline Image getImage(int index) {return ImPyramid[index];};
};

void OpticalFlow::help() {

    printf("-opticalflow computes optical flow between two images based on \n"
           "ECCV 2004 paper from Brox et al.(see source code for credits), \n"
           "which is performed in multiple scales to handle big displacements.\n"
           "This operation returns a three-channel image with x and y pixel\n"
           "offset vectors, and energy cost which can be used as confidence \n"
           "measure for the estimation of each vector. \n"
           "If you have an initial estimate of flow vectors, you give three\n"
           "input images and give any value as an argument to indicate you want to\n"
           "use your initial estimate.\n"
           "\n"
           " arguments [useInitialGuess] \n"
           "  - useInitialGuess (default: none)\n"
           "\n"
           "Usage: ImageStack -load target.jpg -load source.jpg -load guess.tmp -opticalflow 1 -save flow.tmp\n"
           "Usage: ImageStack -load target.jpg -load source.jpg -opticalflow -save flow.tmp\n\n");
}

void OpticalFlow::parse(vector<string> args) {

    if (args.size() > 0) {
        Image result;
        result = apply(stack(1), stack(2), stack(0));
        pop();
        push(result);
    } else {
        Image result;
        result = apply(stack(0), stack(1));
        push(result);
    }
}

Image OpticalFlow::apply(Window source, Window target) {

    // default parameters
    double alpha = 0.005;
    double ratio = 0.75;
    int minWidth = 8;

    int nOuterFPIterations = 20;
    int nInnerFPIterations = 1;
    int nCGIterations = 50;

    Image vx(source.width, source.height, 1, 1);
    Image vy(source.width, source.height, 1, 1);
    Image warpI2(source.width, source.height, source.frames, source.channels);

    Coarse2FineFlow(vx, vy, warpI2, source, target, alpha, ratio, minWidth, nOuterFPIterations, nInnerFPIterations, nCGIterations);

    Image flowVector = Adjoin::apply(vx, vy, 'c');
    Image confidence = evalConfidence(source, warpI2, flowVector, alpha, 1.0);
    Image output = Adjoin::apply(flowVector, confidence, 'c');

    return output;
}

Image OpticalFlow::apply(Window source, Window target, Window initial) {

    // default parameters
    double alpha = 0.005;
    double ratio = 0.75;
    int minWidth = 8;

    int nOuterFPIterations = 20;
    int nInnerFPIterations = 1;
    int nCGIterations = 50;

    printf("Use initial estimate... \n");
    Image temp = Transpose::apply(initial, 'c', 't');
    Window vxInit(temp, 0, 0, 0, initial.width, initial.height, 1);
    Window vyInit(temp, 0, 0, 1, initial.width, initial.height, 1);

    Image vx = Image(vxInit);
    Image vy = Image(vyInit);
    Image warpI2(source.width, source.height, source.frames, source.channels);

    RefineFlow(initial, vx, vy, warpI2, source, target, alpha, ratio, minWidth, nOuterFPIterations, nInnerFPIterations, nCGIterations);

    Image flowVector = Adjoin::apply(vx, vy, 'c');
    Image confidence = evalConfidence(source, warpI2, flowVector, alpha, 1.0);
    Image output = Adjoin::apply(flowVector, confidence, 'c');

    return output;
}

// function to perform coarse to fine optical flow estimation
// notice that vx, vy, warpI2 are output arguments as well as input arguments
void OpticalFlow::Coarse2FineFlow(Image &vx, Image &vy, Image &warpI2, const Image Im1, const Image Im2, double alpha, double ratio, int minWidth, int nOuterFPIterations, int nInnerFPIterations, int nCGIterations) {

    // first build the pyramid of the two images
    GaussianPyramid GPyramid1;
    GaussianPyramid GPyramid2;
    GPyramid1.ConstructPyramid(Im1, ratio, minWidth);
    GPyramid2.ConstructPyramid(Im2, ratio, minWidth);

    // now iterate from the top level to the bottom
    Image Image1, Image2, WarpImage2;

    for (int k = GPyramid1.nlevels()-1; k >= 0; k--) {

        printf("Pyramid level %d ...\n", k);

        int width = GPyramid1.getImage(k).width;
        int height = GPyramid1.getImage(k).height;
        Image1 = im2feature(GPyramid1.getImage(k));
        Image2 = im2feature(GPyramid2.getImage(k));

        if (k == GPyramid1.nlevels()-1) { // if at the top level
            vx = Image(width, height, 1, 1);
            vy = Image(width, height, 1, 1);
            WarpImage2 = Image2.copy();
        } else {
            vx = Resample::apply(vx, width, height);
            Scale::apply(vx, 1.0/ratio);
            vy = Resample::apply(vy, width, height);
            Scale::apply(vy, 1.0/ratio);
            warpFL(WarpImage2,Image1,Image2,vx,vy);
        }

        SmoothFlowPDE(Image1, Image2, WarpImage2, vx, vy, alpha, nOuterFPIterations, nInnerFPIterations, nCGIterations);

    }
    warpFL(warpI2,Im1,Im2,vx,vy);

}

// function to refine given optical flow estimation
// notice that vx, vy, warpI2 are output arguments as well as input arguments
void OpticalFlow::RefineFlow(Image initial, Image &vx, Image &vy, Image &warpI2, const Image Im1, const Image Im2, double alpha, double ratio, int minWidth, int nOuterFPIterations, int nInnerFPIterations, int nCGIterations) {

    Image Image1,Image2,WarpImage2;

    printf("Pyramid level 0 ...\n");

    Image1 = im2feature(Im1);
    Image2 = im2feature(Im2);

    warpFL(WarpImage2,Image1,Image2,vx,vy);

    SmoothFlowPDE(Image1,Image2,WarpImage2,vx,vy,alpha,nOuterFPIterations,nInnerFPIterations,nCGIterations);

    warpFL(warpI2,Im1,Im2,vx,vy);
}

// function to convert image to feature image
Image OpticalFlow::im2feature(Image im) {
    Image imfeature;
    int nchannels=im.channels;

    if (nchannels==1) {

        Image imdx,imdy;

        Image filter(5,1,1,1);
        filter(0,0)[0] = -1;
        filter(1,0)[0] = 8;
        filter(2,0)[0] = 0;
        filter(3,0)[0] = -8;
        filter(4,0)[0] = 1;
        Scale::apply(filter, 1.0/12.0);

        imdx = Convolve::apply(im, filter, Convolve::Clamp, Multiply::Outer);
        imdy = Convolve::apply(im, Transpose::apply(filter, 'x', 'y'), Convolve::Clamp, Multiply::Outer);

        imfeature = Adjoin::apply(im, imdx, 'c');
        imfeature = Adjoin::apply(imfeature, imdy, 'c');

    } else if (nchannels==3) {

        vector<float> grayMatrix;
        for (int i = 0; i < im.channels; i++) {
            grayMatrix.push_back(1.0f/im.channels);
        }
        Image grayImage = ColorMatrix::apply(im, grayMatrix);

        Image imdx,imdy;

        Image filter(5,1,1,1);
        filter(0,0)[0] = -1;
        filter(1,0)[0] = 8;
        filter(2,0)[0] = 0;
        filter(3,0)[0] = -8;
        filter(4,0)[0] = 1;
        Scale::apply(filter, 1.0/12.0);

        imdx = Convolve::apply(grayImage, filter, Convolve::Clamp, Multiply::Outer);
        imdy = Convolve::apply(grayImage, Transpose::apply(filter, 'x', 'y'), Convolve::Clamp, Multiply::Outer);

        // imfeature: gray, dx, dy, G-R, G-B (5 channels)
        imfeature = Image(im.width, im.height, 1, 5);
        float *pGray = grayImage(0,0,0);
        float *pDx = imdx(0,0,0);
        float *pDy = imdy(0,0,0);
        float *pIm = im(0,0,0);
        float *pData = imfeature(0,0,0);

        for (int i=0; i<im.width*im.height; i++) {
            *pData++ = *pGray++;
            *pData++ = *pDx++;
            *pData++ = *pDy++;
            *pData++ = pIm[1] - pIm[0];
            *pData++ = pIm[1] - pIm[2];
            pIm += 3;
        }
    } else {
        imfeature = im.copy();
    }

    return imfeature;
}

// function to compute optical flow field using two fixed point iterations
// Input arguments:
//     Im1, Im2:        frame 1 and frame 2
//  warpIm2:        the warped frame 2 according to the current flow field u and v
//  u,v:            the current flow field, NOTICE that they are also output arguments
void OpticalFlow::SmoothFlowPDE(const Image Im1, const Image Im2, Image &warpIm2, Image &u, Image &v, double alpha, int nOuterFPIterations, int nInnerFPIterations, int nCGIterations) {

    int imWidth,imHeight,nChannels,nPixels;
    imWidth=Im1.width;
    imHeight=Im1.height;
    nChannels=Im1.channels;
    nPixels=imWidth*imHeight;

    Image du(imWidth,imHeight,1,1),dv(imWidth,imHeight,1,1);
    Image uu(imWidth,imHeight,1,1),vv(imWidth,imHeight,1,1);
    Image ux(imWidth,imHeight,1,1),uy(imWidth,imHeight,1,1);
    Image vx(imWidth,imHeight,1,1),vy(imWidth,imHeight,1,1);
    Image Phi_1st(imWidth,imHeight,1,1);
    Image Psi_1st(imWidth,imHeight,1,nChannels);

    Image imdxy,imdx2,imdy2,imdtdx,imdtdy;
    Image ImDxy,ImDx2,ImDy2,ImDtDx,ImDtDy;
    Image A11,A12,A22,b1,b2;
    Image foo1,foo2;

    // variables for conjugate gradient
    Image r1,r2,p1,p2,q1,q2;
    double *rou;
    rou=new double[nCGIterations];

    double varepsilon_phi=pow(0.001,2);
    double varepsilon_psi=pow(0.001,2);

    //--------------------------------------------------------------------------
    // the outer fixed point iteration
    //--------------------------------------------------------------------------
    for (int count=0; count<nOuterFPIterations; count++) {

        Image mask,imdx,imdy,imdt;

        // compute the gradient
        getDxs(imdx,imdy,imdt,Im1,warpIm2);

        // generate the mask to set the weight of the pxiels moving outside of the image boundary to be zero
        genInImageMask(mask,vx,vy);

        // set the derivative of the flow field to be zero
        reset(du);
        reset(dv);

        //--------------------------------------------------------------------------
        // the inner fixed point iteration
        //--------------------------------------------------------------------------
        for (int hh=0; hh<nInnerFPIterations; hh++) {

            // compute the derivatives of the current flow field
            if (hh==0) {
                uu = u.copy();
                vv = v.copy();
            } else {
                uu = add(u,du);
                vv = add(v,dv);
            }

            ux = gradient(uu, 'x');
            uy = gradient(uu, 'y');
            vx = gradient(vv, 'x');
            vy = gradient(vv, 'y');

            // compute the weight of phi
            reset(Phi_1st);
            float *phiData=Phi_1st(0,0,0);
            double temp;
            const float *uxData,*uyData,*vxData,*vyData;
            uxData=ux(0,0,0);
            uyData=uy(0,0,0);
            vxData=vx(0,0,0);
            vyData=vy(0,0,0);
            for (int i=0; i<nPixels; i++) {
                temp=uxData[i]*uxData[i]+uyData[i]*uyData[i]+vxData[i]*vxData[i]+vyData[i]*vyData[i];
                phiData[i]=1.0/(2.0*sqrt(temp+varepsilon_phi));
            }

            // compute the nonlinear term of psi
            reset(Psi_1st);
            float *psiData=Psi_1st(0,0,0);
            const float *imdxData,*imdyData,*imdtData;
            const float *duData,*dvData;
            imdxData=imdx(0,0,0);
            imdyData=imdy(0,0,0);
            imdtData=imdt(0,0,0);
            duData=du(0,0,0);
            dvData=dv(0,0,0);

            //double _a  = 10000, _b = 0.1;
            if (nChannels==1) {
                for (int i=0; i<nPixels; i++) {
                    temp=imdtData[i]+imdxData[i]*duData[i]+imdyData[i]*dvData[i];
                    //if(temp*temp<0.04)
                    psiData[i]=1.0/(2.0*sqrt(temp*temp+varepsilon_psi));
                    //psiData[i] = _a*_b/(1+_a*temp*temp);
                }
            } else {
                for (int i=0; i<nPixels; i++) {
                    for (int k=0; k<nChannels; k++) {
                        int offset=i*nChannels+k;
                        temp=imdtData[offset]+imdxData[offset]*duData[i]+imdyData[offset]*dvData[i];
                        //if(temp*temp<0.04)
                        psiData[offset]=1.0/(2.0*sqrt(temp*temp+varepsilon_psi));
                        //psiData[offset] =  _a*_b/(1+_a*temp*temp);
                    }
                }
            }

            // prepare the components of the large linear system

            ImDxy = multiply(Psi_1st, imdx, imdy);
            ImDx2 = multiply(Psi_1st, imdx, imdx);
            ImDy2 = multiply(Psi_1st, imdy, imdy);
            ImDtDx = multiply(Psi_1st, imdt, imdx);
            ImDtDy = multiply(Psi_1st, imdt, imdy);

            if (nChannels>1) {
                imdxy = collapse(ImDxy);
                imdx2 = collapse(ImDx2);
                imdy2 = collapse(ImDy2);
                imdtdx = collapse(ImDtDx);
                imdtdy = collapse(ImDtDy);
            } else {
                imdxy = ImDxy.copy();
                imdx2 = ImDx2.copy();
                imdy2 = ImDy2.copy();
                imdtdx = ImDtDx.copy();
                imdtdy = ImDtDy.copy();
            }

            // filtering
            Image filter(3,1,1,1);
            filter(0,0)[0] = 0.2;
            filter(1,0)[0] = 0.6;
            filter(2,0)[0] = 0.2;

            A11 = Convolve::apply(imdx2, filter, Convolve::Clamp, Multiply::Outer);
            A11 = Convolve::apply(A11, Transpose::apply(filter, 'x', 'y'), Convolve::Clamp, Multiply::Outer);
            A12 = Convolve::apply(imdxy, filter, Convolve::Clamp, Multiply::Outer);
            A12 = Convolve::apply(A12, Transpose::apply(filter, 'x', 'y'), Convolve::Clamp, Multiply::Outer);
            A22 = Convolve::apply(imdy2, filter, Convolve::Clamp, Multiply::Outer);
            A22 = Convolve::apply(A22, Transpose::apply(filter, 'x', 'y'), Convolve::Clamp, Multiply::Outer);

            // add epsilon to A11 and A22
            Offset::apply(A11, alpha*0.1);
            Offset::apply(A22, alpha*0.1);

            // form b
            b1 = Convolve::apply(imdtdx, filter, Convolve::Clamp, Multiply::Outer);
            b1 = Convolve::apply(b1, Transpose::apply(filter, 'x', 'y'), Convolve::Clamp, Multiply::Outer);
            b2 = Convolve::apply(imdtdy, filter, Convolve::Clamp, Multiply::Outer);
            b2 = Convolve::apply(b2, Transpose::apply(filter, 'x', 'y'), Convolve::Clamp, Multiply::Outer);

            // laplacian filtering of the current flow field
            Laplacian(foo1,u,Phi_1st);
            Laplacian(foo2,v,Phi_1st);

            float *b1Data,*b2Data;
            const float *foo1Data,*foo2Data;
            b1Data=b1(0,0,0);
            b2Data=b2(0,0,0);
            foo1Data=foo1(0,0,0);
            foo2Data=foo2(0,0,0);

            for (int i=0; i<nPixels; i++) {
                b1Data[i] = -b1Data[i] - alpha*foo1Data[i];
                b2Data[i] = -b2Data[i] - alpha*foo2Data[i];
            }

            //-----------------------------------------------------------------------
            // conjugate gradient algorithm
            //-----------------------------------------------------------------------
            r1 = b1.copy();
            r2 = b2.copy();
            reset(du);
            reset(dv);

            for (int k=0; k<nCGIterations; k++) {

                rou[k]=norm2(r1)+norm2(r2);
                //cout<<rou[k]<<endl;
                if (rou[k]<1E-10) {
                    break;
                }

                if (k==0) {
                    p1 = r1.copy();
                    p2 = r2.copy();
                } else {
                    double ratio=rou[k]/rou[k-1];
                    p1 = addAfterScale(r1, p1, ratio);
                    p2 = addAfterScale(r2, p2, ratio);
                }

                // go through the large linear system
                foo1 = multiply(A11, p1);
                foo2 = multiply(A12, p2);
                q1 = add(foo1, foo2);
                Laplacian(foo1,p1,Phi_1st);
                q1 = addAfterScale(q1, foo1, alpha);

                foo1 = multiply(A12, p1);
                foo2 = multiply(A22, p2);
                q2 = add(foo1, foo2);
                Laplacian(foo2,p2,Phi_1st);
                q2 = addAfterScale(q2, foo2, alpha);

                double beta;
                beta=rou[k]/(innerproduct(p1,q1)+innerproduct(p2,q2));

                du = addAfterScale(du, p1, beta);
                dv = addAfterScale(dv, p2, beta);

                r1 = addAfterScale(r1, q1, -beta);
                r2 = addAfterScale(r2, q2, -beta);

            }
            //-----------------------------------------------------------------------
            // end of conjugate gradient algorithm
            //-----------------------------------------------------------------------
            //printf("end cg\n");


        }// end of inner fixed point iteration

        // the following procedure is merely for debugging
        //cout<<"du "<<norm2(du)<<" dv "<<norm2(dv)<<endl;

        // update the flow field
        Add::apply(u, du);
        Add::apply(v, dv);



        warpFL(warpIm2,Im1,Im2,u,v);

    }// end of outer fixed point iteration

    delete rou;
}

void OpticalFlow::warpFL(Image &warpIm2, Image Im1, Image Im2, Image vx, Image vy) {

    // resample Im2 based on vx, vy offsets
    warpIm2 = Image(Im2.width, Im2.height, 1, Im2.channels);
    float *pWarp = warpIm2(0,0,0);
    float *pVx = vx(0,0,0);
    float *pVy = vy(0,0,0);

    for (int y=0; y<Im2.height; y++) {
        for (int x=0; x<Im2.width; x++) {
            Im2.sample2D(x + *pVx++, y + *pVy++, pWarp);
            pWarp += Im2.channels;
        }
    }
}

//  function to compute dx, dy and dt for motion estimation
void OpticalFlow::getDxs(Image &imdx, Image &imdy, Image &imdt, Image im1, Image im2) {

    // Im1 and Im2 are the smoothed version of im1 and im2
    Image Im1 = im1.copy(), Im2 = im2.copy();
    FastBlur::apply(Im1, 1.0, 1.0, 0);
    FastBlur::apply(Im2, 1.0, 1.0, 0);

    Image filter(5,1,1,1);
    filter(0,0)[0] = -1;
    filter(1,0)[0] = 8;
    filter(2,0)[0] = 0;
    filter(3,0)[0] = -8;
    filter(4,0)[0] = 1;
    Scale::apply(filter, 1.0/12.0);

    imdx = Convolve::apply(Im2, filter, Convolve::Clamp, Multiply::Outer);
    imdy = Convolve::apply(Im2, Transpose::apply(filter, 'x', 'y'), Convolve::Clamp, Multiply::Outer);
    imdt = subtract(Im2, Im1);
}

// function to generate mask of the pixels that move inside the image boundary
void OpticalFlow::genInImageMask(Image &mask, Image vx, Image vy) {

    int imWidth,imHeight;
    imWidth=vx.width;
    imHeight=vx.height;

    if (mask.width!=vx.width || mask.height!=vx.height || mask.frames!=vx.frames || mask.channels!=vx.channels) {
        mask = Image(vx.width, vx.height, vx.frames, vx.channels);
    }

    float *pVx = vx(0,0,0);
    float *pVy = vy(0,0,0);
    float *pMask = mask(0,0,0);
    double x,y;
    for (int i=0; i<imHeight; i++) {
        for (int j=0; j<imWidth; j++) {
            y = i + *pVx++;
            x = j + *pVy++;
            if (x<0  || x>imWidth-1 || y<0 || y>imHeight-1) {
                *pMask++ = 0;
            } else {
                *pMask++ = 1;
            }
        }
    }
}

void OpticalFlow::reset(Image &im) {
    float *ptr = im(0,0,0);
    for (int i=0; i<im.width*im.height*im.frames*im.channels; i++) {
        *ptr++ = 0;
    }
}

Image OpticalFlow::collapse(Image in) {
    Image out(in.width, in.height, 1, 1);
    float *inPtr = in(0,0,0);
    float *outPtr = out(0,0,0);
    float scale = 1.0 / in.channels;
    for (int i=0; i<in.width*in.height; i++) {
        float temp = 0;
        for (int j=0; j<in.channels; j++) {
            temp += *inPtr++;
        }
        *outPtr++ = temp * scale;
    }
    return out;
}

void OpticalFlow::Laplacian(Image &output, Image input, Image weight) {

    output = Image(input.width, input.height, 1, 1);
    reset(output);

    const float *inputData=input(0,0,0),*weightData=weight(0,0,0);
    int width=input.width,height=input.height;
    Image foo(width,height,1,1);
    float *fooData=foo(0,0,0),*outputData=output(0,0,0);

    // horizontal filtering
    for (int i=0; i<height; i++) {
        for (int j=0; j<width-1; j++) {
            int offset=i*width+j;
            fooData[offset]=(inputData[offset+1]-inputData[offset])*weightData[offset];
        }
    }
    for (int i=0; i<height; i++) {
        for (int j=0; j<width; j++) {
            int offset=i*width+j;
            if (j<width-1) {
                outputData[offset]-=fooData[offset];
            }
            if (j>0) {
                outputData[offset]+=fooData[offset-1];
            }
        }
    }
    reset(foo);

    // vertical filtering
    for (int i=0; i<height-1; i++) {
        for (int j=0; j<width; j++) {
            int offset=i*width+j;
            fooData[offset]=(inputData[offset+width]-inputData[offset])*weightData[offset];
        }
    }
    for (int i=0; i<height; i++) {
        for (int j=0; j<width; j++) {
            int offset=i*width+j;
            if (i<height-1) {
                outputData[offset]-=fooData[offset];
            }
            if (i>0) {
                outputData[offset]+=fooData[offset-width];
            }
        }
    }
}

double OpticalFlow::norm2(Image im) {

    float *ptr = im(0,0,0);
    double temp = 0;

    for (int i=0; i<im.width*im.height*im.frames*im.channels; i++) {
        temp += ((*ptr) * (*ptr));
        ptr++;
    }

    return temp;
}

double OpticalFlow::innerproduct(Image im1, Image im2) {

    float *ptr1 = im1(0,0,0);
    float *ptr2 = im2(0,0,0);
    double temp = 0;

    for (int i=0; i<im1.width*im1.height*im1.frames*im1.channels; i++) {
        temp += *ptr1++ * *ptr2++;
    }

    return temp;
}

// return a + b
Image OpticalFlow::add(Image a, Image b) {

    Image output = Image(a.width, a.height, a.frames, a.channels);
    float *pOutput = output(0,0,0);
    float *pA = a(0,0,0);
    float *pB = b(0,0,0);

    for (int i=0; i<a.width*a.height*a.frames*a.channels; i++) {
        *pOutput++ = *pA++ + *pB++;
    }

    return output;
}

// return a - b
Image OpticalFlow::subtract(Image a, Image b) {

    Image output = Image(a.width, a.height, a.frames, a.channels);
    float *pOutput = output(0,0,0);
    float *pA = a(0,0,0);
    float *pB = b(0,0,0);

    for (int i=0; i<a.width*a.height*a.frames*a.channels; i++) {
        *pOutput++ = *pA++ - *pB++;
    }

    return output;
}

Image OpticalFlow::gradient(Image im, char dimension) {

    Image output(im.width, im.height, 1, im.channels);
    int minx = 0, miny = 0;
    int dx = 0, dy = 0;

    if (dimension == 'x') {
        dx = 1; minx = 1;
    } else if (dimension == 'y') {
        dy = 1; miny = 1;
    }

    // walk backwards through the data, looking at the untouched data for the differences
    for (int y = im.height - 1; y >= miny; y--) {
        for (int x = im.width - 1; x >= minx; x--) {
            for (int c = 0; c < im.channels; c++) {
                output(x, y, 0)[c] = im(x, y, 0)[c] - im(x - dx, y - dy, 0)[c];
            }
        }
    }
    return output;
}

// return a * b
Image OpticalFlow::multiply(Image a, Image b) {

    Image output = Image(a.width, a.height, a.frames, a.channels);
    float *pOutput = output(0,0,0);
    float *pA = a(0,0,0);
    float *pB = b(0,0,0);

    for (int i=0; i<a.width*a.height*a.frames*a.channels; i++) {
        *pOutput++ = *pA++ * *pB++;
    }

    return output;
}

// return a * b * c
Image OpticalFlow::multiply(Image a, Image b, Image c) {

    Image output = Image(a.width, a.height, a.frames, a.channels);
    float *pOutput = output(0,0,0);
    float *pA = a(0,0,0);
    float *pB = b(0,0,0);
    float *pC = c(0,0,0);

    for (int i=0; i<a.width*a.height*a.frames*a.channels; i++) {
        *pOutput++ = *pA++ * *pB++ * *pC++;
    }

    return output;
}

// return a + b * c
Image OpticalFlow::addAfterScale(Image a, Image b, double c) {

    Image output = Image(a.width, a.height, a.frames, a.channels);
    float *pOutput = output(0,0,0);
    float *pA = a(0,0,0);
    float *pB = b(0,0,0);

    for (int i=0; i<a.width*a.height*a.frames*a.channels; i++) {
        *pOutput++ = *pA++ + *pB++ * c;
    }

    return output;
}

Image OpticalFlow::phi(Image a) {

    Image output = Image(a.width, a.height, a.frames, a.channels);
    float *pOutput = output(0,0,0);
    float *pA = a(0,0,0);
    float eps = 0.001*0.001;

    for (int i=0; i<a.width*a.height*a.frames*a.channels; i++) {
        *pOutput++ = sqrt(*pA * *pA + eps);
        pA++;
    }

    return output;
}

Image OpticalFlow::evalConfidence(Image source, Image target, Image flow, float alpha, float gamma) {

    Image data(source.width, source.height, 1, 1);
    Image smooth(source.width, source.height, 1, 1);
    Image sourceDx = gradient(source, 'x');
    Image sourceDy = gradient(source, 'y');
    Image targetDx = gradient(target, 'x');
    Image targetDy = gradient(target, 'y');
    Image flowDx = gradient(flow, 'x');
    Image flowDy = gradient(flow, 'y');

    float *Iw = new float[source.channels];
    float *IwDx = new float[source.channels];
    float *IwDy = new float[source.channels];
    for (int y=0; y<source.height; y++) {
        for (int x=0; x<source.width; x++) {
            target.sample2D(x, y, Iw);
            targetDx.sample2D(x, y, IwDx);
            targetDy.sample2D(x, y, IwDy);
            for (int c=0; c<source.channels; c++) {
                data(x, y)[0] += (Iw[c] - source(x, y)[c]) * (Iw[c] - source(x, y)[c]);
                data(x, y)[0] += gamma * (IwDx[c] - sourceDx(x, y)[c]) * (IwDx[c] - sourceDx(x, y)[c]);
                data(x, y)[0] += gamma * (IwDy[c] - sourceDy(x, y)[c]) * (IwDy[c] - sourceDy(x, y)[c]);
            }
            smooth(x, y)[0] += (flowDx(x, y)[0] * flowDx(x, y)[0]);
            smooth(x, y)[0] += (flowDx(x, y)[1] * flowDx(x, y)[1]);
            smooth(x, y)[0] += (flowDy(x, y)[0] * flowDy(x, y)[0]);
            smooth(x, y)[0] += (flowDy(x, y)[1] * flowDy(x, y)[1]);
        }
    }
    delete Iw; delete IwDx; delete IwDy;

    Image confidence = addAfterScale(phi(data), phi(smooth), alpha);

    return confidence;
}

void OpticalFlowWarp::help() {

    printf("-opticalflowwarp warps an image to match the other image based on \n"
           "optical flow estimation. \n"
           "\n"
           "Usage: ImageStack -load from.jpg -load to.jpg -opticalflowwarp -save warped.jpg\n\n"
          );
}

void OpticalFlowWarp::parse(vector<string> args) {

    assert(args.size() == 0, "-opticalflowwarp takes no arguments\n");

    Image result;
    result = apply(stack(1), stack(1), stack(0));
    push(result);
}

Image OpticalFlowWarp::apply(Window input, Window from, Window to) {

    Image flow = OpticalFlow::apply(to, from);

    vector<string> str;
    str.push_back(string("(x+[0])/width"));
    str.push_back(string("(y+[1])/height"));
    Image coord = EvalChannels::apply(flow, str);

    Image output = Warp::apply(coord, input);

    return output;
}

#include "footer.h"

