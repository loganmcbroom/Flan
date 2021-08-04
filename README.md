# Flan
A non-real-time audio processing library based on the phase vocoder transform.

# Building

### Dependencies:

#### libsndfile:
```
git clone https://github.com/libsndfile/libsndfile.git
cd libsndfile
cmake -S . -B build -DCMAKE_DEBUG_POSTFIX=d -DCMAKE_BUILD_TYPE=Release
cmake --build build --target install --config Release
```
There is also a windows installer available here: http://www.mega-nerd.com/libsndfile/#Download


#### FFTWF:
See https://www.fftw.org/

Building FFTW is OS-dependant. Unix systems can download the source from the FFTW website (http://www.fftw.org/fftw-3.3.9.tar.gz) and run:
```
./configure
make
make install
```

Windows systems can follow the instructions here (http://www.fftw.org/install/windows.html) to install an import library. 
The FFTW devs seem to really dislike Windows, but the instructions work.

Either way, the official github repo should be avoided, as you need to run a preprocessor on it to generate the actual source code.

#### OpenCL:
This should exist on modern computers already.

### Flan:
```
git clone https://github.com/loganmcbroom/flan
cd flan
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target install --config Release
```
The "-DCMAKE_BUILD_TYPE=Release" is for single-config generators, while "--config Release" is for multi-config generators.
Flan also exports the build tree, so "--target install" isn't needed if the build will remain where it was built.

# Usage
The classes Audio and PVOC represent all the main algorithms in flan. They inherit buffer functionality from AudioBuffer and PVOCBuffer. 
The Function class and its children allow a great deal of freedom in passing functions (usually in the form of lambdas) in place of constants 
when calling Audio and PVOC methods. Some basic Audio synthesis functions are found in Synthesis.h.

Any questions on installation or usage can be directed to:

Discord: https://discord.gg/QmQuJKB

Email: loganmcbroom@gmail.com
