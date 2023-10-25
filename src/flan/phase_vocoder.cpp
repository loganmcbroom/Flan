#include "phase_vocoder.h"

namespace flan {

MF phase_vocoder( double & phase_buffer, std::complex<float> cpx, Frequency bin_frequency, FrameRate analysis_rate, FrameRate sample_rate )
	{
	/*
	The phase vocoder algorithm can be applied to any time-frequency transform with complex output that uses time-overlapping analysis windows.
	It looks at how the phase of the complex outputs of that transform changes over time within a given output bin. By comparing that change
	to the change in phase that is expected for the bin, the frequency of the partial that is being detected can be predicted.
	
	Details on phase wrapping:
		When using an analysis rate slower than the audio sample rate, there is an extra step called "phase wrapping". Phase wrapping forces
		a particular angle called "delta phase" to be on the interval [-pi,pi]. The computations below show the possible frequency outputs
		for both using and not using phase wrapping.
		
		Without phase wrapping:
			frequency
			= bin_frequency + delta_frequency 
			= bin_frequency + delta_phase * analysis_rate / pi2 )
			= bin_frequency + ( phase_diff - expected_phase_diff ) * analysis_rate / pi2
			= bin_frequency + ( phase_diff - bin_frequency / analysis_rate * pi2 ) * analysis_rate / pi2
			= bin_frequency + phase_diff * analysis_rate / pi2 ) - bin_frequency
			= phase_diff * analysis_rate / pi2 )
			= [-pi2, pi2] * analysis_rate / pi2
			= [-analysis_rate, analysis_rate] <- fine for analysis rate = sample rate, otherwise this becomes a low-pass filter

		With phase wrapping:
			frequency
			= bin_frequency + delta_frequency 
			= bin_frequency + wrapped_delta_phase * analysis_rate / pi2
			= bin_frequency + [-pi,pi] * analysis_rate / pi2
			= bin_frequency + [-analysis_rate/2,analysis_rate/2] <- Allows full frequency range even for small analysis rates.

		Phase wrapping, then, must happen if the analysis rate is less than the sample rate.
	*/
	const bool use_wrapping = analysis_rate < sample_rate;
	auto wrap = []( Radian x ) -> Radian
		{
		return x - pi2 * std::round( x / pi2 );
		};

	const Radian phase = std::arg( cpx );
	const Radian phase_diff = phase - phase_buffer; // Change in phase from the last analysis time
	phase_buffer = phase; // Update phase buffer

	const Radian expected_phase_diff = bin_frequency / analysis_rate * pi2; // Expected phase change
	const Radian delta_phase = phase_diff - expected_phase_diff; // The difference between the actual and expected phase change
	const Radian wrapped_delta_phase = use_wrapping ? wrap( delta_phase ) : delta_phase; // Wrap delta phase into [-pi,pi]
	const Frequency delta_frequency = wrapped_delta_phase * analysis_rate / pi2;

	return { std::abs( cpx ), bin_frequency + delta_frequency };
	}

std::complex<float> inverse_phase_vocoder( double & phase_buffer, MF mf, FrameRate analysis_rate )
	{
	const Radian phase_diff = mf.f / analysis_rate * pi2;
	phase_buffer += phase_diff;
	if( phase_buffer > pi2 ) phase_buffer = std::fmod( phase_buffer, pi2 );
	return std::polar( mf.m, float( phase_buffer ) );
	}
	
}