Ideas:
	Falter:
		only allow premade interpolators, lua ones aren't thread safe
			Consider not ever taking interpolators

	continuity maintenance in wavelength/frequency finders
	octave error minimization in yin-based pitch finder (fast contour finder)
	logarithmic specral analysis https://essentia.upf.edu/reference/std_LogSpectrum.html
	convolve constant signals via fft
	contour pitch mean should be magnitude weighted
	denormals?: https://www.carlh.net/plugins/denormals.php
	ir generation - https://dsp.stackexchange.com/questions/41696/calculating-the-inverse-filter-for-the-exponential-sine-sweep-method/41700#41700
	Audio/PV conversions can compute all non-overlapping windows at the same time using this http://fftw.org/fftw3_doc/Advanced-Complex-DFTs.html
	perturb sucks currently but might have potential using something like simplex noise, see comments on process
		it also could be placing frequencies outside bins that can hold them, use the bin shift strategy as in PV::extrapolateTime
	PV::synthesize should be using this: https://zynaddsubfx.sourceforge.io/doc_0.html
		there are a bunch of other places this can be applied
		prism, where it says "// There is no harmonic to scale, so create one"
	PV::timeExtrapolate could use a Func2x1 instead of interp, not sure the best way to set it up though

Process ideas:
	Audio::removeSilence
	Audio Dynamics
	Audio::splitChannels, Audio::joinChannels
	audio to function of frequency
	Audio::select using PSOLA
	PV::alignHarmonics
	PV::shape variant that passes time, bin frequency, and the MF to the shaper
	PV::compress for some type of multiband compressor
	Graph::convertToPV
	Drum sequencing tools
	Wavetable::PADsynth https://zynaddsubfx.sourceforge.io/doc/PADsynth/PADsynth.htm
	Paul Stretch - this is just stft stretching with random phase
		This algorithm may be improvable

Task:
	fix buffer_access name conflict (free function vs method)
		just use the PV name for the other ones

	find a better way to handle function sampling, then remove buffer_access

	compressor:
		need to use some kind of stereo linking, is there any reason to use anything besides max?
			For now, no.
		Harmonic distortion needs to be minimal while creating an acceptable imitation of instantaneous control
		look-ahead is a non issue with a proper nrt design
		The idea that you need a fast attack is just a consequence or rt systems
			It may just be cultural perception, but I think people expect fast attack and slow exponential decay in dB domain


	
	PV::select is misplacing frequencies, is there a fix for this?
		yes, using whatever strategy to get a frequency from the selected location, multiply by destBinFreq/selectedFreq
		just rewrite select
		test just using nearest neighbor MF instead of max-mag

	Some PV algorithms may be worth porting to other formats, either that or they could be rewritten to work for all those formats
		PV is fast but requires MF combining
		SPV is slow 
		SQPV has mid speed, but seems to have some kind of time stretch issue? Check that out.
			lol

	SPV and SQPV graphing, PV too actually, needs to handle that bins don't hold the frequency at that bin all the time

	Pitch is still worth exploring, we can just have Pitch be an object that converts to and from frequency
		This would be a lot nicer if not everything could implicit convert to frequency



	filter package
		1 pole
			-multimode
			-low
			-high
			-low shelf
			-high shelf
		SFV
			-multimode
			-low 
			-band
			-high
			notch
			band shelf and peak have a chapter 10 section
		arbitrary order butterworth
			-lowpass
			-highpass

		arbitrary order band shelf using chapter 10. The 1-pole shelves are base cases of this and can be removed once this is done.
		
		Ladder?
		Parks-McLellen design
		Linkwitz-Riley for splitting into bands
		check out non-linear filters
