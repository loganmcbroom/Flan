# flan
Creative audio processing library.

# Compilation
Optionally depends on libsndfile, FFTW (the float version), and OpenCL.
libsndfile enables file formats besides WAVE.
FFTW makes several main algorithms significantly faster.
OpenCL permits large optimizations in some of the main algorithms.
Things will hopefully work out of the box via cmake. 
Dependencies are cmake options and enabled by default, disabling them should allow building anywhere c++14 can be used.
The cmake find modules are fairly bare-bones, so if they can't find something make sure it is in the system PATH.

# Usage
The classes Audio and PVOC represent all the main algorithms in flan. They inherit buffer functionality from AudioBuffer and PVOCBuffer. 
The Function class and its children allow a great deal of freedom in passing functions (usually in the form of lambdas) in place of constants 
when calling Audio and PVOC methods. Some basic Audio synthesis functions are found in Synthesis.h. 
Most other usages can be found in the documentation, here: https://codedocs.xyz/loganmcbroom/flan/

Any questions on installation or usage can be directed to:

Discord: https://discord.gg/QmQuJKB

Email: loganmcbroom@gmail.com
