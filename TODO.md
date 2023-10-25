Ideas:
	continuity maintenance in wavelength/frequency finders
	octave error minimization in yin-based pitch finder (fast contour finder)
	contour pitch mean should be magnitude weighted
	ir generation - https://dsp.stackexchange.com/questions/41696/calculating-the-inverse-filter-for-the-exponential-sine-sweep-method/41700#41700
	Audio/PV conversions can compute all non-overlapping windows at the same time using this http://fftw.org/fftw3_doc/Advanced-Complex-DFTs.html
	perturb sucks currently but might have potential using something like simplex noise, see comments on process
		it also could be placing frequencies outside bins that can hold them, use the bin shift strategy as in PV::extrapolateTime
	Some procs should be using this: https://zynaddsubfx.sourceforge.io/doc_0.html
		prism where it says "// There is no harmonic to scale, so create one"
	PV::time_extrapolate could use a Func2x1 instead of interp, not sure the best way to set it up though
	PV::select is misplacing frequencies, is there a fix for this?
		yes, using whatever strategy to get a frequency from the selected location, multiply by destBinFreq/selectedFreq
		just rewrite select
		test just using nearest neighbor MF instead of max-mag
	SPV and SQPV graphing, PV too actually, needs to handle that bins don't hold the frequency at that bin all the time
	find a better way to handle function sampling, then remove buffer_access
	fix perturb

Process ideas:
	Audio to function of frequency
	PV::alignHarmonics
	PV::shape variant that passes time, bin frequency, and the MF to the shaper
	PV::compress for some type of multiband compressor
	Graph::convert_to_PV
	Drum sequencing tools
	Wavetable::PADsynth https://zynaddsubfx.sourceforge.io/doc/PADsynth/PADsynth.htm
	Paul Stretch - this is just stft stretching with random phase
		This algorithm may be improvable
	Audio::limiter
		look-ahead option in compressor would take care of this

Task:
	Audio::split_at_frequencies

	disable AudioMod copyable