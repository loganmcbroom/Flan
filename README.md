# XCDP
Creative audio processing library based on the Composers Desktop Project.

Link against libsndfile, libsamplerate, and fftw. Fftw requires a shared library, at least on windows where I live.

I'm new to cmake so this isn't exported yet (TODO)

# Usage:
	There are only two classes you need to be aware of, Audio and PVOC. They inherit buffer functionality from AudioBuffer and PVOCBuffer. 
	Pretty simple. Construct an Audio with a file name and that file is loaded into the Audio. Audio methods process their object and return the result.
	Note that all methods are const. If you want to work in the spectral domain call convertToPVOC, convertToAudio gets you back.
	Any questions on installation or usage can be directed to:
	Discord: https://discord.gg/QmQuJKB
	Email: loganmcbroom@gmail.com
