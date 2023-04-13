#include "phase_vocoder.h"

#include "Utility/MF.h"

namespace flan {

MF phase_vocoder( Radian & phase_buffer, std::complex<float> cpx, Frequency bin_frequency, Time hop_time )
	{
	const Radian phase = std::arg( cpx );

	const Radian phase_diff = phase - phase_buffer; // Actual change in radians for the current hop
	phase_buffer = phase; // Update phase buffer

	const Radian expected_phase_diff = bin_frequency * hop_time * pi2; // Expected change in radians
	const Radian delta_phase = phase_diff - expected_phase_diff; // The difference between the actual and the expected phase changes, aka how far off was our sinusoid from the expected.
	const Radian wrapped_delta_phase = delta_phase - pi2 * std::round( delta_phase / pi2 ); // Wrap delta phase into [-pi,pi]
	const Frequency delta_frequency = wrapped_delta_phase / ( hop_time * pi2 );

	return { std::abs( cpx ), bin_frequency + delta_frequency };
	}

std::complex<float> inverse_phase_vocoder( Radian & phase_buffer, MF mf, Time hop_time )
	{
	const Radian phase_diff = mf.f * hop_time * pi2;
	phase_buffer += phase_diff;
	return std::polar( mf.m, phase_buffer );
	}
	
}