toilet
kitchen cabs
mirrors
microwave
palish appliances/clean air fryer
trash maintain
bathroom counter/sink
tub
sweep/vac
mop
put away dishes maintain

Maybe todo:
	contour pitch mean should be magnitude weighted
	ir generation - https://dsp.stackexchange.com/questions/41696/calculating-the-inverse-filter-for-the-exponential-sine-sweep-method/41700#41700
	Audio/PV conversions can compute all non-overlapping windows at the same time using this http://fftw.org/fftw3_doc/Advanced-Complex-DFTs.html
	perturb sucks currently but might have potential using something like simplex noise, see comments on process
		it also could be placing frequencies outside bins that can hold them, use the bin shift strategy as in PV::time_extrapolate
	PV::time_extrapolate could use a Func2x1 instead of interp, not sure the best way to set it up though

Process ideas:
	Audio::limiter
		look-ahead option in compressor would take care of this
	Audio::split_at_frequencies for multiple splits
	PV::align_harmonics
	PV::chorus might work, based on expanding partials
		Per frame we should be able to detect partials and expand them
		The main issue I see is partial continuity 

Task:
	Something's in the Lake at Bagel Woods

	in prism, this: "bound_bin( frequency_to_bin( freq ) - 10 );" and the following line use 10 arbitraryly
		This seems bad to use without testing, but it might be fine
		notes_close may also be a bad metric

	factor get_pointers to... somewhere

	ohhhh maaaaaan time_to_frame always rounding down might be causing extremely small issues all over

	explodes:
		synth.filter_1pole_multinotch( 5, []( Second t ){ return 10000*t; }, .9 );


