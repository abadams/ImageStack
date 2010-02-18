#!/usr/bin/bash
wget http://savannah.nongnu.org/download/openexr/openexr-1.4.0a.tar.gz
wget http://www.fftw.org/fftw-3.1.2.tar.gz
wget http://www.libsdl.org/release/SDL-1.2.11.tar.gz

for DEP in *.tar.gz; do
    tar xvzf $DEP
done

rm *.tar.gz

echo "To install OpenEXR and fftw, just do the usual:"
echo "./configure && make && make install"
echo 
echo "When configuring SDL on win32, you should probably use"
echo "./configure --disable-stdio-redirect"
echo "to prevent it from stealing stdin and stderr."
