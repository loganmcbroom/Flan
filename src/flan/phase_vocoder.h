#pragma once

#include <complex>

#include "defines.h"
#include "flan/Utility/MF.h"

namespace flan {

// The phase buffer should be a Radian type, but it needs additional precision in this case to avoid rounding error accumulation.

MF phase_vocoder( double & phase_buffer, std::complex<float> cpx, Frequency bin_frequency, FrameRate analysis_rate, FrameRate sample_rate );
std::complex<float> inverse_phase_vocoder( double & phase_buffer, MF mf, FrameRate analysis_rate );

}