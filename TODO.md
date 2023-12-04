Maybe todo:
	contour pitch mean should be magnitude weighted
	ir generation - https://dsp.stackexchange.com/questions/41696/calculating-the-inverse-filter-for-the-exponential-sine-sweep-method/41700#41700
	Audio/PV conversions can compute all non-overlapping windows at the same time using this http://fftw.org/fftw3_doc/Advanced-Complex-DFTs.html
	perturb sucks currently but might have potential using something like simplex noise, see comments on process
		it also could be placing frequencies outside bins that can hold them, use the bin shift strategy as in PV::extrapolateTime
	PV::time_extrapolate could use a Func2x1 instead of interp, not sure the best way to set it up though
	find a better way to handle function sampling, then remove buffer_access
	Speed up stereo spatialize
	find a general method to avoid sampling constant functions

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
	some way to load midi file as Function

	"Instrument" concept
		synthesizes Audio from notes
		midi analog format is a set of notes
		Like wavetables, this is an object that is first built up, then used for synthesis
		Really it's just a converter from Note to Audio
		the converter function can decide how to output audio based on the note
		the basic utility is to simply return some fixed audio pitched and scaled as in a sampler
		This is sounding a lot like synth_grains with no jitter
		The main change is that a note has some definitive parameters
		I don't want to be overly general and lose out on specific interesting ideas like each adsr parameter being controllable via note params
			That should be relatively simple