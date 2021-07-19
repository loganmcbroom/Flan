# Flan
Creative audio processing library.

# Building
Dependancies:
	libsndfile:
		git clone https://github.com/libsndfile/libsndfile.git
		cd libsndfile
		cmake -S . -B build
		cmake --build build --target install
		There is also a windows installer available here: http://www.mega-nerd.com/libsndfile/#Download

	FFTWF:
		https://www.fftw.org/
		I wish you luck.

	OpenCL: 
		This should exist on modern computers already.

Flan:
	git clone https://github.com/loganmcbroom/flan
	cd flan
	cmake -S . -B build
	cmake --build build --target install --config Release

# Usage
The classes Audio and PVOC represent all the main algorithms in flan. They inherit buffer functionality from AudioBuffer and PVOCBuffer. 
The Function class and its children allow a great deal of freedom in passing functions (usually in the form of lambdas) in place of constants 
when calling Audio and PVOC methods. Some basic Audio synthesis functions are found in Synthesis.h. 
Most other usages can be found in the documentation, here: https://codedocs.xyz/loganmcbroom/flan/

Codedocs is having trouble updating at the moment.

Any questions on installation or usage can be directed to:

Discord: https://discord.gg/QmQuJKB

Email: loganmcbroom@gmail.com
