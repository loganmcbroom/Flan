Ideas:
	Falter:
		only allow premade interpolators, lua ones aren't thread safe

	continuity maintenance in wavelength/frequency finders
	octave error minimization in yin-based pitch finder (fast contour finder)
	logarithmic specral analysis https://essentia.upf.edu/reference/std_LogSpectrum.html
	convolve constant signals via fft
	contour pitch mean should be magnitude weighted
	denormals?: https://www.carlh.net/plugins/denormals.php
	ir generation - https://dsp.stackexchange.com/questions/41696/calculating-the-inverse-filter-for-the-exponential-sine-sweep-method/41700#41700
	Audio/PVOC conversions can compute all non-overlapping windows at the same time using this http://fftw.org/fftw3_doc/Advanced-Complex-DFTs.html
	perturb sucks currently but might have potential using something like simplex noise, see comments on process
		it also could be placing frequencies outside bins that can hold them, use the bin shift strategy as in PVOC::extrapolateTime
	pvoc::synthesize should be using this: https://zynaddsubfx.sourceforge.io/doc_0.html
		there are a bunch of other places this can be applied
		prism, where it says "// There is no harmonic to scale, so create one"
	PVOC::timeExtrapolate could use a Func2x1 instead of interp, not sure the best way to set it up though

Process ideas:
	pvoc::alignHarmonics
	Graph::convertToPVOC
	PVOC::shape variant that passes time, bin frequency, and the MF to the shaper
	Audio::removeSilence
	Drum sequencing tools
	Audio Dynamics
	Wavetable::PADsynth https://zynaddsubfx.sourceforge.io/doc/PADsynth/PADsynth.htm
	Audio::splitChannels, Audio::joinChannels
	audio to function of frequency

Task:
	filter package

	compressor, doesn't need to be perfect but at least do some research

	PSOLA package
		It sounds cool, but what are the goals?
			Coherent stretching
			Potentially a data type, as in Audio::convertToPSOLA which would chop the Audio into grains with some kind of precomputed overlap data.
				If so, what procs?
					Stretch - this is a special case of select
					Select
		If I can't think of anything besides stretch and select, these should just be Audio algorithms
	
	select is misplacing frequencies, is there a fix for this?
		yes, using whatever strategy to get a frequency from the selected location, multiply by destBinFreq/selectedFreq
		just rewrite select
		test just using nearest neighbor MF instead of max-mag

	PV
		rename PVOC to PV?
	SPV
	SQPV
		
	
			