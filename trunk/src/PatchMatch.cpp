
// CS448F final project
// Implementation of PatchMatch algorithm and its applications
// Sung Hee Park (shpark7@stanford.edu)


#include "main.h"
#include "File.h"
#include "Geometry.h"
#include "PatchMatch.h"
#include "Arithmetic.h"
#include "Calculus.h"
#include "Statistics.h"
#include "Filter.h"

// PATCHMATCH =============================================================//

void PatchMatch::help() {

    printf("-patchmatch computes approximate nearest neighbor field from the top\n"
	   "image on the stack to the second image on the stack, using the\n"
	   "algorithm from the PatchMatch SIGGRAPH 2009 paper. This operation\n"
	   "requires two input images which may have multiple frames.\n"
	   "It returns an image with four channels. First three channels \n"
	   "correspond to x, y, t coordinate of closest patch and \n"
	   "fourth channels contains the sum of squared differences \n"
	   "between patches. \n"
	   "\n"
	   " arguments [numIter] [patchSize]\n"
	   "  - numIter : number of iterations performed. (default: 5)\n"
	   "  - patchSize : size of patch. (default: 7, 7x7 square patch)\n"
	   " You can omit some arguments from right to use default values.\n"
	   "\n"
	   "Usage: ImageStack -load target.jpg -load source.jpg -patchmatch -save match.tmp\n\n");
}

void PatchMatch::parse(vector<string> args) {

    int numIter = 5, patchSize = 7;

    assert(args.size() <= 2, "-patchmatch takes two or fewer arguments\n");
    if (args.size() == 2) {
	numIter = readInt(args[0]);
	patchSize = (int) readInt(args[1]);
    } else if (args.size() == 1) {
	numIter = readInt(args[0]);
    }

    Image result;		     

    result = apply(stack(0), stack(1), numIter, patchSize);

    push(result);
}

Image PatchMatch::apply(Window source, Window target, int iterations, int patchSize) {
    return apply(source, target, Window(), iterations, patchSize);
}

Image PatchMatch::apply(Window source, Window target, Window mask, int iterations, int patchSize) {

    if (mask) {
	assert(target.width == mask.width &&
	       target.height == mask.height &&
	       target.frames == mask.frames, 
	       "Mask must have the same dimensions as the target\n");
	assert(mask.channels == 1,
	       "Mask must have a single channels\n");
    }
    assert(iterations > 0, "Iterations must be a strictly positive integer\n");
    assert(patchSize >= 3 && (patchSize & 1), "Patch size must be at least 3 and odd\n");

    // convert patch diameter to patch radius
    patchSize /= 2;

    // For each source pixel, output a 3-vector to the best match in
    // the target, with an error as the last channel.
    Image out(source.width, source.height, source.frames, 4);
    
    // Iterate over source frames, finding a match in the target where
    // the mask is high
    
    float *outPtr = out(0, 0, 0);
    for (int t = 0; t < source.frames; t++) {	       
	// INITIALIZATION - uniform random assignment
	for(int y=0; y < source.height; y++) {
	    for(int x=0; x < source.width; x++) {
		int dx = rand() % target.width;
		int dy = rand() % target.height;
		int dt = rand() % target.frames;
		*outPtr++ = rand() % target.width;
		*outPtr++ = rand() % target.height;
		*outPtr++ = rand() % target.frames;
		*outPtr++ = distance(source, target, mask,
				     t, x, y,
				     dx, dy, dt,
				     patchSize, HUGE_VAL);
	    }
	}
    }

    bool forwardSearch = true;

    for (int i=0; i < iterations; i++) {

	//printf("Iteration %d\n", i);
	
	// PROPAGATION
	if(forwardSearch) {
	    
	    // Forward propagation - compare left, center and up
	    for (int t = 0; t < source.frames; t++) {
		for(int y = 1; y < source.height; y++) {
		    outPtr = out(1, y, t);
		    float *leftPtr = out(0, y, t);
		    float *upPtr = out(0, y-1, t);
		    for(int x = 1; x < source.width; x++) {
			float distLeft = distance(source, target, mask, 
						  t, x, y, 
						  leftPtr[2], leftPtr[0]+1, leftPtr[1], 
						  patchSize, outPtr[3]);
			
			if (distLeft < outPtr[3]) {
			    outPtr[0] = leftPtr[0]+1;
			    outPtr[1] = leftPtr[1];
			    outPtr[2] = leftPtr[2];
			    outPtr[3] = distLeft;
			}

			float distUp = distance(source, target, mask, 
						t, x, y, 
						upPtr[2], upPtr[0], upPtr[1]+1, 
						patchSize, outPtr[3]);

			if (distUp < outPtr[3]) {
			    outPtr[0] = upPtr[0];
			    outPtr[1] = upPtr[1]+1;
			    outPtr[2] = upPtr[2];
			    outPtr[3] = distUp;			    
			}

			outPtr += 4;
			leftPtr += 4;
			upPtr += 4;

			// TODO: Consider searching across time as well
		    }
		}
	    }

	} else {		
	    
	    // Backward propagation - compare right, center and down
	    for (int t = source.frames-1; t >= 0; t--) {
		for(int y = source.height-2; y >= 0; y--) {
		    outPtr = out(source.width-2, y, t);
		    float *rightPtr = out(source.width-1, y, t);
		    float *downPtr = out(source.width-2, y+1, t);
		    for(int x = source.width-2; x >= 0; x--) {
			float distRight = distance(source, target, mask,
						   t, x, y, 
						   rightPtr[2], rightPtr[0]-1, rightPtr[1],
						   patchSize, outPtr[3]);
			
			if (distRight < outPtr[3]) {
			    outPtr[0] = rightPtr[0]-1;
			    outPtr[1] = rightPtr[1];
			    outPtr[2] = rightPtr[2];
			    outPtr[3] = distRight;
			}

			float distDown = distance(source, target, mask,
						  t, x, y, 
						  downPtr[2], downPtr[0], downPtr[1]-1, 
						  patchSize, outPtr[3]);

			if (distDown < outPtr[3]) {
			    outPtr[0] = downPtr[0];
			    outPtr[1] = downPtr[1]-1;
			    outPtr[2] = downPtr[2];
			    outPtr[3] = distDown;			    
			}

			outPtr -= 4;
			rightPtr -= 4;
			downPtr -= 4;
			// TODO: Consider searching across time as well
			
		    }
		}
	    }
	}	    

	forwardSearch = !forwardSearch;

	// RANDOM SEARCH
	float *outPtr = out(0, 0, 0);

	for (int t = 0; t < source.frames; t++) {
	    for(int y = 0; y < source.height; y++) {
		for(int x = 0; x < source.width; x++) {

		    int radius = target.width > target.height ? target.width : target.height;
		    
		    // search an exponentially smaller window each iteration
		    while (radius > 1) {
			// Search around current offset vector (distance-weighted)
			
			// clamp the search window to the image
			int minX = (int)outPtr[0] - radius;
			int maxX = (int)outPtr[0] + radius + 1;
			int minY = (int)outPtr[1] - radius;
			int maxY = (int)outPtr[1] + radius + 1;
			if (minX < 0) minX = 0;
			if (maxX > target.width) maxX = target.width;
			if (minY < 0) minY = 0;
			if (maxY > target.height) maxY = target.height;
			
			int randX = rand() % (maxX - minX) + minX;
			int randY = rand() % (maxY - minY) + minY;
			int randT = rand() % target.frames;
			float dist = distance(source, target, mask,
					      t, x, y,
					      randT, randX, randY,
					      patchSize, outPtr[3]);
			if (dist < outPtr[3]) {
			    outPtr[0] = randX;
			    outPtr[1] = randY;
			    outPtr[2] = randT;
			    outPtr[3] = dist;
			}
			
			radius >>= 1;
			
		    }		    
		    outPtr += 4;
		}
	    }
	}
    }

    return out;
}

float PatchMatch::distance(Window source, Window target, Window mask,
			   int st, int sx, int sy, 
			   int tt, int tx, int ty,
			   int patchSize, float prevDist) {

    // Do not use patches on boundaries
    if (tx < patchSize || tx >= target.width-patchSize || 
	ty < patchSize || ty >= target.height-patchSize) {
	return HUGE_VAL;
    }

    // Compute distance between patches
    // Average L2 distance in RGB space
    float dist = 0;
    int count = 0;

    float threshold = prevDist*target.channels*(2*patchSize+1)*(2*patchSize+1);

    int x1 = max(-patchSize, -sx, -tx);
    int x2 = min(patchSize, -sx+source.width-1, -tx+target.width-1);
    int y1 = max(-patchSize, -sy, -ty);
    int y2 = min(patchSize, -sy+source.height-1, -ty+target.height-1);
    
    /*
    int x1 = -patchSize, x2 = patchSize;
    int y1 = -patchSize, y2 = patchSize;
    */

    for(int y = y1; y <= y2; y++) {

	float *pSource = source(sx+x1, sy+y, st);
	float *pTarget = target(tx+x1, ty+y, tt);
	float *pMask = NULL;
	if (mask) pMask = mask(tx+x1, ty+y, tt);

	for (int i = 0; i <= x2-x1; i++) {
	    float d = 0;
	    for (int j = 0; j < target.channels; j++) {
		d += (*pSource - *pTarget)*(*pSource - *pTarget);
		count++;
		pSource++; pTarget++;
	    }

	    // Divide by the mask weight here
	    
	    if (mask) {
		if (*pMask == 0) return HUGE_VAL;
		d /= (*pMask);
		pMask += mask.channels;
	    }

	    dist += d;

	    // Early termination
	    if (dist > threshold) {return HUGE_VAL;}
	}
    }

    if (!count) return HUGE_VAL;

    return dist / count;
}


// BIDIRECTIONAL SIMILARITY =====================================================//


void BidirectionalSimilarity::help() {
    printf("-bidirectional minimizes bidirectional similarity measure \n"
	   "based on patchmatch algorithm for acceleration. This is operated\n"
	   "only on the finest scale. The target image is assumed to have\n"
	   "one frame for now.\n"
	   "\n"
	   " arguments [alpha] [numIter] [numIterPM]\n"
	   "  - alpha : weighting between completeness and coherence term. (default: 0.5)\n"
	   "  - numIter : number of iterations performed. (default: 5)\n"
	   "  - numIterPM : number of iterations performed in PatchMatch. (default: 5)\n"
	   " You can omit some arguments from right to use default values.\n"
	   "\n" 
	   "Usage: ImageStack -load source.jpg -load target.jpg -bidirectional 0.5 -display\n\n");
}

void BidirectionalSimilarity::parse(vector<string> args) {

    float alpha = 0.5;
    int numIter = 5;
    int numIterPM = 5;

    assert(args.size() <= 3, "-bidirectional takes three or fewer arguments\n");
    if (args.size() == 3) {
 	alpha = readFloat(args[0]);
	numIter = readFloat(args[1]);
	numIterPM = readFloat(args[2]);
    } else if (args.size() == 2) {
	alpha = readFloat(args[0]);
	numIter = readFloat(args[1]);
    } else if (args.size() == 1) {
	alpha = readFloat(args[0]);
    } 
    
    apply(stack(1), stack(0), Window(), Window(), alpha, numIter, numIterPM);
}


// Reconstruct the portion of the target where the mask is high, using
// the portion of the source where its mask is high. Source and target
// masks are allowed to be null windows.
void BidirectionalSimilarity::apply(Window source, Window target, 
				    Window sourceMask, Window targetMask,
				    float alpha, int numIter, int numIterPM) {

    

    // TODO: intelligently crop the input to where the mask is high +
    // patch radius on each side


    // recurse
    if (source.width > 32 && source.height > 32 && target.width > 32 && target.height > 32) {
	Image smallSource = Resample::apply(source, source.width/2, source.height/2, source.frames);
	Image smallTarget = Resample::apply(target, target.width/2, target.height/2, target.frames);
	
	Image smallSourceMask;
	Image smallTargetMask;
	if (sourceMask) {
	    smallSourceMask = Resample::apply(sourceMask, 
					      sourceMask.width/2, 
					      sourceMask.height/2,
					      sourceMask.frames);
	}

	if (targetMask) {
	    smallTargetMask = Resample::apply(targetMask,
					      targetMask.width/2, 
					      targetMask.height/2, 
					      targetMask.frames);
	}

	apply(smallSource, smallTarget, smallSourceMask, smallTargetMask, alpha, numIter, numIterPM);

	Image newTarget = Resample::apply(smallTarget, target.width, target.height, target.frames);
	
	for (int t = 0; t < target.frames; t++) {
	    for (int y = 0; y < target.height; y++) {
		float *targPtr = target(0, y, t);
		float *newTargPtr = newTarget(0, y, t);
		memcpy(targPtr, newTargPtr, sizeof(float)*target.channels*target.width);
	    }
	}	
    }

    for(int i = 0; i < numIter; i++) {
	printf("Bidir.. %d/%d\n", i, numIter);

	int patchSize = 5; 
	Image completeMatch, coherentMatch;

	// The homogeneous output for this iteration
	Image out(target.width, target.height, target.frames, target.channels+1);

	printf("Completeness...\n");

 	if (alpha != 0) {

	    // COMPLETENESS TERM
	    Image completeMatch = PatchMatch::apply(source, target, targetMask, numIterPM, patchSize);

	    // For every patch in the source, splat it onto the
	    // nearest match in the target, weighted by the source
	    // mask and also by the inverse of the patch distance
	    float *matchPtr = completeMatch(0, 0, 0);
	    for (int t = 0; t < source.frames; t++) {
		for (int y = 0; y < source.height; y++) {
		    float *srcMaskPtr = sourceMask(0, y, t);
		    for (int x = 0; x < source.width; x++) {

			if (!sourceMask || srcMaskPtr[0] > 0) {
			
			    int dstX = (int)matchPtr[0];
			    int dstY = (int)matchPtr[1];
			    int dstT = (int)matchPtr[2];			
			    float weight = 1.0f/(matchPtr[3]+1);			
			    
			    if (sourceMask) weight *= srcMaskPtr[0];
			    
			    for (int dy = -patchSize/2; dy <= patchSize/2; dy++) {
				if (y+dy < 0) continue;
				if (y+dy >= source.height) break;
				float *sourcePtr = source(x-patchSize/2, y+dy, t);
				float *outPtr = out(dstX-patchSize/2, dstY+dy, dstT);
				for (int dx = -patchSize/2; dx <= patchSize/2; dx++) {
				    if (x+dx < 0) {
					outPtr += out.channels;
					sourcePtr += source.channels;
				    } else if (x+dx >= source.width) {
					break;
				    } else {
					for (int c = 0; c < source.channels; c++) {
					    (*outPtr++) += weight*(*sourcePtr++);
					}
					(*outPtr++) += weight;
				    }
				}
			    }
			}

			srcMaskPtr++;
			matchPtr += completeMatch.channels;
		    }
		}
	    }
	}
	
	printf("Coherence...\n"); fflush(stdout);

	if (alpha != 1) {

	    // COHERENCE TERM
	    Image coherentMatch = PatchMatch::apply(target, source, sourceMask,
						    numIterPM, patchSize);

	    // For every patch in the target, pull from the nearest match in the source
	    float *matchPtr = coherentMatch(0, 0, 0);
	    for (int t = 0; t < target.frames; t++) {
		for (int y = 0; y < target.height; y++) {
		    float *targMaskPtr = targetMask(0, y, t);
		    for (int x = 0; x < target.width; x++) {

			if (!targetMask || targMaskPtr[0] > 0) {
			
			    int dstX = (int)matchPtr[0];
			    int dstY = (int)matchPtr[1];
			    int dstT = (int)matchPtr[2];			
			    float weight = 1.0f/(matchPtr[3]+1);			
			    
			    if (targetMask) weight *= targMaskPtr[0];
			    
			    for (int dy = -patchSize/2; dy <= patchSize/2; dy++) {
				if (y+dy < 0) continue;
				if (y+dy >= out.height) break;
				float *sourcePtr = source(dstX-patchSize/2, dstY+dy, dstT);
				float *outPtr = out(x-patchSize/2, y+dy, t);
				for (int dx = -patchSize/2; dx <= patchSize/2; dx++) {
				    if (x+dx < 0) {
					outPtr += out.channels;
					sourcePtr += source.channels;
				    } else if (x+dx >= out.width) {
					break;
				    } else {
					for (int c = 0; c < source.channels; c++) {
					    (*outPtr++) += weight*(*sourcePtr++);
					}
					(*outPtr++) += weight;
				    }
				}
			    }			    
			}

			targMaskPtr++;
			matchPtr += coherentMatch.channels;
		    }
		}
	    }	   
	}
	
	printf("Reconstructing...\n"); fflush(stdout);

	// rewrite the target using the homogeneous output
	float *outPtr = out(0, 0, 0);
	float *targMaskPtr = targetMask(0, 0, 0);
	for (int t = 0; t < out.frames; t++) {
	    for (int y = 0; y < out.height; y++) {
		float *targetPtr = target(0, y, t);
		for (int x = 0; x < out.width; x++) {
		    float w = 1.0f/(outPtr[target.channels]);
		    float a = 1;
		    if (targetMask) {
			 a = *targMaskPtr++;
		    }
		    if (a == 1) {
			for (int c = 0; c < target.channels; c++) {
			    targetPtr[0] = w*(*outPtr++);
			    targetPtr++;
			}
		    } else if (a > 0) {
			for (int c = 0; c < target.channels; c++) {
			    targetPtr[0] *= (1-a);
			    targetPtr[0] += a*w*(*outPtr++);
			    targetPtr++;
			}
		    } 
		    outPtr++;
		}
	    }
	}
	
    }
}
