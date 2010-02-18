#define NO_MAIN

// includes that don't survive well in the namespace
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <list>
#include <set>

#ifdef WIN32
#include <winsock2.h>
#else
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <netdb.h>
#endif

#ifndef NO_OPENEXR
#include <ImfRgbaFile.h>
#include <ImfArray.h>
#endif

#include <math.h>

namespace ImageStack {

#define IMAGESTACK_NAMESPACE

#include "main.cpp"
#include "Calculus.cpp"
#include "Color.cpp"
#include "Control.cpp"
#include "Convolve.cpp"
#include "Complex.cpp"
#include "Display.cpp"
#include "DisplayWindow.cpp"
#include "DFT.cpp"
#include "Exception.cpp"
#include "File.cpp"
#include "FileCSV.cpp"
#include "FileEXR.cpp"
#include "FileFLO.cpp"
#include "FileHDR.cpp"
#include "FileJPG.cpp"
#include "FilePNG.cpp"
#include "FilePPM.cpp"
#include "FileTGA.cpp"
#include "FileTIFF.cpp"
#include "FileTMP.cpp"
#include "FileWAV.cpp"
#include "Filter.cpp"
#include "Geometry.cpp"
#include "HDR.cpp"
#include "Image.cpp"
#include "LAHBPCG.cpp"
#include "LightField.cpp"
#include "Arithmetic.cpp"
#include "Network.cpp"
#include "NetworkOps.cpp"
#include "Operation.cpp"
#include "Paint.cpp"
#include "Panorama.cpp"
#include "PatchMatch.cpp"
#include "Parser.cpp"
#include "Prediction.cpp"
#include "Projection.cpp"
#include "Stack.cpp"
#include "Statistics.cpp"
#include "Wavelet.cpp"
#include "Alignment.cpp"
#include "GaussTransform.cpp"
#include "WLS.cpp"

}


