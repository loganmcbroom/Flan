#pragma once

#include <complex>
#include "defines.h"

namespace flan {

struct MF;

// This executes the crux of a phase vocoder transform, that is estimating the frequency 
// of a partial using phase information over time
MF phase_vocoder( Radian & phase_buffer, std::complex<float> cpx, Frequency bin_frequency, Time hop_time );

std::complex<float> inverse_phase_vocoder( Radian & phase_buffer, MF mf, Time hop_time );

}