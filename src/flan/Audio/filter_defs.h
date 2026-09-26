#pragma once

#include <complex>
#include "flan/Function.h"

using namespace flan;
using Pole = std::complex<float>;
using Mix_1pole = std::array<float,2>;
using Mix_2pole = std::array<float,3>;
using Mix_Func_1pole = Function<Second, Mix_1pole>;
using Mix_Func_2pole = Function<Second, Mix_2pole>;

struct Filter_1Pole 
    {
	Filter_1Pole( FrameRate sr );
	std::array<Sample, 2> process_sample( Sample x, Frequency cutoff_unwarped, bool use_prewarp = true );
	Sample process_sample_and_mix( Sample x, Frequency cutoff_unwarped, Mix_1pole mix, bool use_prewarp = true );

	Sample s;
	const float T_half;
    };

struct Filter_2Pole {
	Filter_2Pole( FrameRate sr );
	std::array<Sample, 3> process_sample( Sample x, Frequency cutoff_unwarped, float R, bool use_prewarp = true );
	Sample process_sample_and_mix( Sample x, Frequency cutoff_unwarped, float R, Mix_2pole mix, bool use_prewarp = true );

	Sample s1, s2;
	const float T_half;
};