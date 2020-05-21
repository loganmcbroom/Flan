# XCDP
Creative audio processing library based on the Composers Desktop Project.

# Compilation
Depends on libsndfile, libsamplerate, fftw (the float version). 
OpenCL is an optional dependancy that permits large optimizations in some of the main algorithms.
Things will hopefully work out of the box via cmake.
The cmake find modules are fairly bare-bones, so if they can't find something make sure it is in the system PATH.

# Usage
The classes Audio and PVOC represent all the main algorithms in XCDP. They inherit buffer functionality from AudioBuffer and PVOCBuffer. 
The Function class and its children allow a great deal of freedom in passing functions (usually in the form of lambdas) in place of constants 
when calling Audio and PVOC methods. Some basic Audio synthesis functions are found in Synthesis.h. 
Most other usages can be found in the documentation, here: https://codedocs.xyz/loganmcbroom/XCDP/

Any questions on installation or usage can be directed to:

Discord: https://discord.gg/QmQuJKB

Email: loganmcbroom@gmail.com
