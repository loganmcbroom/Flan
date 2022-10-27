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
There is also a windows installer available [here](http://www.mega-nerd.com/libsndfile/#Download).


#### FFTWF:
See https://www.fftw.org/

Building FFTW is OS-dependant. Unix systems can download the source from the [FFTW website](http://www.fftw.org/fftw-3.3.9.tar.gz) and run:
```
./configure
make
make install
```

Windows systems can follow the instructions [here](http://www.fftw.org/install/windows.html) to install an import library. 
The FFTW devs seem to really dislike Windows, but the instructions work.

Either way, the official github repo should be avoided, as you need to run a preprocessor on it to generate the actual source code.

#### OpenCL:
This should exist on most modern computers already. If OpenCL is missing from your system, there are a number of distributions of it from a number of vendors.
Installing any of them should supply your system with a FindOpenCL.cmake file which can find all the needed files. Flan uses c++ bindings that not all distributions contain.
If that is the case with your distribution, download the binding repository [here](https://github.com/KhronosGroup/OpenCL-CLHPP), and copy the opencl.hpp and CL2.hpp files to the same folder as opencl.h.

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
The classes flan::Audio and flan::PVOC represent all the main algorithms in flan. They inherit buffer functionality from flan::AudioBuffer and flan::PVOCBuffer. The flan::Function class and its children allow a great deal of freedom in passing functions (usually in the form of lambdas) in place of constants when calling flan::Audio and flan::PVOC methods. Some basic synthesis functions are found in flan::Synthesis.

Any further questions can be directed to:

Discord: https://discord.gg/QmQuJKB

Email: loganmcbroom@gmail.com
