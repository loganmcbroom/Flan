Ideas:
	filter package
	continuity maintenance in wavelength/frequency finders
	octave error minimization in yin-based pitch finder (fast contour finder)
	audio to function of frequency
	convolve constant signals via fft
	Audio::removeSilence
	Audio::cut uses an extra buffer copy when calling fades, remove it
	contour pitch mean should be magnitude weighted
	only generate docs for flan stuff
	pvoc::alignHarmonics
	select/mix/join should consider bitrate conversion when needed
	PVOC::shape variant that passes time, bin frequency, and the MF to the shaper
	can PVOC::freeze be made with a more normal input type?

Bugs:
	compute_d array access crash in AudioInformation - bounds check
	select seems slow, check it out

Task:
	Expand main page documentation - maybe examples