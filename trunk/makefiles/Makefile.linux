BIN_CCFLAGS = -std=gnu++0x -O3 -Winvalid-pch -Wshadow -Wall -Werror -Wno-uninitialized -pipe -march=native -ffast-math -fopenmp -rdynamic


LIB_CCFLAGS = $(BIN_CCFLAGS) -fPIC

##########################
# Dependencies           #
##########################

# Comment out any dependencies you do not wish to include, or do not have
SDL_CCFLAGS = -I /usr/include/SDL
SDL_LIBS = `sdl-config --libs`

JPEG_CCFLAGS =
JPEG_LIBS = -ljpeg

TIFF_CCFLAGS =
TIFF_LIBS = -ltiff

PNG_CCFLAGS = 
PNG_LIBS = -lpng

FFTW_CCFLAGS = 
FFTW_LIBS = -lfftw3f

#OPENEXR_CCFLAGS = -I /usr/local/include/OpenEXR
#OPENEXR_LIBS = -L/usr/local/lib -lImath -lHalf -lIex -lIlmImf

#################################
# Target names and installation #
#################################

BIN_TARGET = bin/ImageStack
default: $(BIN_TARGET)
LIB_TARGET = lib/libImageStack.so
LIB_FLAGS = -shared -fPIC

install: $(BIN_TARGET)
	cp $(BIN_TARGET) /usr/local/bin
	cp src/imagestack /etc/bash_completion.d/
	bash -e ./src/imagestack

libinstall: library
	cp $(LIB_TARGET) /usr/local/lib
	mkdir -p /usr/local/include/ImageStack
	cp include/*.h /usr/local/include/ImageStack/

include Makefile.common











