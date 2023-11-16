Maybe todo:
	contour pitch mean should be magnitude weighted
	ir generation - https://dsp.stackexchange.com/questions/41696/calculating-the-inverse-filter-for-the-exponential-sine-sweep-method/41700#41700
	Audio/PV conversions can compute all non-overlapping windows at the same time using this http://fftw.org/fftw3_doc/Advanced-Complex-DFTs.html
	perturb sucks currently but might have potential using something like simplex noise, see comments on process
		it also could be placing frequencies outside bins that can hold them, use the bin shift strategy as in PV::extrapolateTime
	PV::time_extrapolate could use a Func2x1 instead of interp, not sure the best way to set it up though
	find a better way to handle function sampling, then remove buffer_access
	Speed up stereo spatialize

Process ideas:
	Drum sequencing tools
	Wavetable::PADsynth https://zynaddsubfx.sourceforge.io/doc/PADsynth/PADsynth.htm
	Paul Stretch - this is just stft stretching with random phase
		This algorithm may be improvable
	Audio::limiter
		look-ahead option in compressor would take care of this
	Audio::split_at_frequencies for multiple splits
	PV::align_harmonics

Task: